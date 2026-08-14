# Random-effects structure builders for BBmm (RI / RI+RS, diag or correlated).

#' Normalize corr argument to canonical "unstructured" or "diag"
#' @keywords internal
.normalize_corr <- function(corr) {
  if (is.null(corr)) return("unstructured")
  if (is.numeric(corr) || is.integer(corr)) {
    return(if (as.integer(corr)[1] == 0L) "diag" else "unstructured")
  }
  c0 <- tolower(as.character(corr)[1])
  # Canonical + aliases (cor/correlated kept as informal shortcuts)
  if (c0 %in% c("unstructured", "us", "un", "full",
                "cor", "correlated")) {
    return("unstructured")
  }
  if (c0 %in% c("diag", "diagonal", "ind", "independent", "uncorrelated")) {
    return("diag")
  }
  stop(
    "corr must be 'unstructured'/'us' or 'diag' ",
    "(aliases: cor, correlated, independent, ...)",
    call. = FALSE
  )
}

#' @keywords internal
.is_corr_diag <- function(corr) {
  identical(.normalize_corr(corr), "diag")
}

#' @keywords internal
.n_theta_sigma <- function(q, corr) {
  q <- as.integer(q)
  if (.is_corr_diag(corr)) return(q)
  q * (q + 1L) / 2L
}

#' @keywords internal
.theta0_sigma <- function(q, corr, sd0 = 0.5) {
  q <- as.integer(q)
  if (.is_corr_diag(corr)) {
    return(rep(log(sd0), q))
  }
  # Chol lower: diag = log(sd), off-diag = 0 -> Sigma = diag(sd^2)
  th <- numeric(q * (q + 1L) / 2L)
  k <- 1L
  for (i in seq_len(q)) {
    for (j in seq_len(i)) {
      th[k] <- if (i == j) log(sd0) else 0
      k <- k + 1L
    }
  }
  th
}

#' Unpack Sigma (and sd, corr) from unconstrained theta
#' @keywords internal
.sigma_from_theta <- function(theta, q, corr) {
  q <- as.integer(q)
  if (.is_corr_diag(corr)) {
    sd <- exp(as.numeric(theta))
    Sigma <- diag(sd^2, nrow = q)
    Corr <- diag(q)
    return(list(Sigma = Sigma, sd = sd, Corr = Corr))
  }
  L <- matrix(0, q, q)
  k <- 1L
  th <- as.numeric(theta)
  for (i in seq_len(q)) {
    for (j in seq_len(i)) {
      L[i, j] <- if (i == j) exp(th[k]) else th[k]
      k <- k + 1L
    }
  }
  Sigma <- L %*% t(L)
  sd <- sqrt(pmax(diag(Sigma), 0))
  Corr <- diag(q)
  for (i in seq_len(q)) {
    for (j in seq_len(q)) {
      Corr[i, j] <- Sigma[i, j] / (sd[i] * sd[j] + 1e-12)
    }
  }
  list(Sigma = Sigma, sd = sd, Corr = Corr)
}

#' Find '|'-terms in a random-effects formula (simplified lme4::findbars)
#' @keywords internal
.findbars <- function(term) {
  if (is.null(term) || is.name(term) || is.numeric(term) || is.character(term)) {
    return(NULL)
  }
  if (!is.call(term) && !is.language(term)) return(NULL)
  op <- term[[1L]]
  if (identical(op, as.name("|"))) return(list(term))
  if (identical(op, as.name("("))) return(.findbars(term[[2L]]))
  if (identical(op, as.name("+"))) {
    return(c(.findbars(term[[2L]]), .findbars(term[[3L]])))
  }
  # ~ rhs
  if (identical(op, as.name("~"))) {
    if (length(term) == 3L) return(.findbars(term[[3L]]))
    return(.findbars(term[[2L]]))
  }
  NULL
}

#' Build Z columns for one (lhs | group) term, interleaved by group
#' @keywords internal
.make_Z_bar <- function(lhs, group_fac, data, term_names = NULL) {
  group_fac <- droplevels(as.factor(group_fac))
  G <- nlevels(group_fac)
  lev <- levels(group_fac)
  mm <- stats::model.matrix(lhs, data = data)
  q <- ncol(mm)
  if (is.null(term_names)) term_names <- colnames(mm)
  n <- nrow(mm)
  Z <- matrix(0, n, G * q)
  cn <- character(G * q)
  for (g in seq_len(G)) {
    rows <- which(group_fac == lev[g])
    for (t in seq_len(q)) {
      col <- (g - 1L) * q + t
      Z[rows, col] <- mm[rows, t]
      cn[col] <- paste0(lev[g], ":", term_names[t])
    }
  }
  colnames(Z) <- cn
  list(Z = Z, G = G, q = q, term_names = term_names, group_levels = lev,
       group_name = NULL)
}

#' Normalize random= / random.formula= / Z= into block structure
#'
#' @return list with Z, blocks (list), theta0, data fields for TMB, names, etc.
#' @keywords internal
.build_re_structure <- function(random = NULL, random.formula = NULL,
                                Z = NULL, nRandComp = NULL,
                                data, corr = "unstructured",
                                nObs) {
  corr_default <- .normalize_corr(corr)
  blocks <- list()
  Z_list <- list()

  # ---- Path 1: explicit Z (legacy / advanced) ----
  if (!is.null(Z)) {
    if (!is.null(random) || !is.null(random.formula)) {
      stop("Specify only one of random / random.formula / Z", call. = FALSE)
    }
    if (is.null(nRandComp)) {
      stop("nRandComp must be specified when Z is given", call. = FALSE)
    }
    Z <- as.matrix(Z)
    nRandComp <- as.integer(nRandComp)
    if (ncol(Z) != sum(nRandComp)) {
      stop("sum(nRandComp) must equal ncol(Z)", call. = FALSE)
    }
    # Interpret as independent q=1 diagonal blocks (PROreg parity),
    # unless a single block with re_terms attribute.
    re_terms <- attr(Z, "re_terms")
    re_corr <- attr(Z, "re_corr")
    if (!is.null(re_terms) && length(nRandComp) == 1L) {
      # One grouping factor, ncol(Z)=G*q, nRandComp = G*q or G?
      q <- as.integer(re_terms)
      if (ncol(Z) %% q != 0L) stop("ncol(Z) not divisible by re_terms", call. = FALSE)
      G <- ncol(Z) / q
      corr_b <- if (is.null(re_corr)) corr_default else .normalize_corr(re_corr)
      if (q == 1L) corr_b <- "diag"
      blocks[[1]] <- list(
        name = "1", G = as.integer(G), q = q, corr = corr_b,
        term_names = paste0("term", seq_len(q))
      )
      Z_list[[1]] <- Z
    } else {
      # Each nRandComp entry = one RI factor (q=1)
      col0 <- 0L
      for (i in seq_along(nRandComp)) {
        G <- nRandComp[i]
        Zi <- Z[, col0 + seq_len(G), drop = FALSE]
        col0 <- col0 + G
        blocks[[i]] <- list(
          name = as.character(i), G = G, q = 1L, corr = "diag",
          term_names = "(Intercept)"
        )
        Z_list[[i]] <- Zi
      }
    }
  } else if (!is.null(random)) {
    # ---- Path 2: random = ~ (1+time|id) or list(id = ~ 1 + time) ----
    if (!is.null(random.formula)) {
      stop("Specify random or random.formula, not both", call. = FALSE)
    }
    if (is.list(random) && !inherits(random, "formula")) {
      # list(id = ~ 1 + time, site = ~ 1)
      nm <- names(random)
      if (is.null(nm) || any(!nzchar(nm))) {
        stop("random list must be named by grouping factors", call. = FALSE)
      }
      for (i in seq_along(random)) {
        gname <- nm[i]
        lhs <- random[[i]]
        if (!inherits(lhs, "formula")) {
          stop("random[['", gname, "']] must be a formula like ~ 1 + time",
               call. = FALSE)
        }
        if (length(lhs) == 3L) {
          stop("Use ~ 1 + time (no LHS response) in random list", call. = FALSE)
        }
        gfac <- eval(as.name(gname), envir = data, enclos = parent.frame())
        bar <- .make_Z_bar(lhs, gfac, data)
        q <- bar$q
        corr_b <- if (q == 1L) "diag" else corr_default
        blocks[[i]] <- list(
          name = gname, G = bar$G, q = q, corr = corr_b,
          term_names = bar$term_names, group_levels = bar$group_levels
        )
        Z_list[[i]] <- bar$Z
      }
    } else if (inherits(random, "formula")) {
      bars <- .findbars(random)
      if (!length(bars)) {
        # No '|': treat as legacy ~ f1 + f2 intercepts
        return(.build_re_structure(
          random.formula = random, data = data, corr = corr_default, nObs = nObs
        ))
      }
      for (i in seq_along(bars)) {
        bcall <- bars[[i]]
        # (lhs | rhs)
        lhs <- bcall[[2L]]
        rhs <- bcall[[3L]]
        # wrap lhs as formula
        if (identical(lhs, 1) || identical(lhs, 1L)) {
          lhs_f <- ~ 1
        } else {
          lhs_f <- eval(parse(text = paste0("~", paste(deparse(lhs), collapse = ""))))
        }
        gname <- all.vars(rhs)
        if (length(gname) != 1L) {
          stop("Grouping factor on the right of '|' must be a single variable",
               call. = FALSE)
        }
        # Fix: when group_fac is evaluated from name in data frame column
        gfac <- tryCatch(
          eval(rhs, envir = as.list(data), enclos = baseenv()),
          error = function(e) data[[gname]]
        )
        if (is.null(gfac)) gfac <- data[[gname]]
        bar <- .make_Z_bar(lhs_f, gfac, data)
        q <- bar$q
        corr_b <- if (q == 1L) "diag" else corr_default
        blocks[[i]] <- list(
          name = gname, G = bar$G, q = q, corr = corr_b,
          term_names = bar$term_names, group_levels = bar$group_levels
        )
        Z_list[[i]] <- bar$Z
      }
    } else {
      stop("random must be a formula or a named list of formulas", call. = FALSE)
    }
  } else if (!is.null(random.formula)) {
    # ---- Path 3: legacy random.formula = ~ z (+ w) ----
    rf <- stats::update(random.formula, ~ . - 1)
    random.mf <- stats::model.frame(formula = rf, data = data)
    nComp <- ncol(random.mf)
    for (i in seq_len(nComp)) {
      gfac <- random.mf[[i]]
      gname <- names(random.mf)[i]
      Zi <- stats::model.matrix(~ gfac - 1)
      colnames(Zi) <- paste0(gname, levels(as.factor(gfac)))
      blocks[[i]] <- list(
        name = gname, G = ncol(Zi), q = 1L, corr = "diag",
        term_names = "(Intercept)",
        group_levels = levels(as.factor(gfac))
      )
      Z_list[[i]] <- Zi
    }
  } else {
    stop("Random part must be specified (random, random.formula, or Z)",
         call. = FALSE)
  }

  Z <- do.call(cbind, Z_list)
  if (nrow(Z) != nObs) stop("nrow(Z) must equal length(y)", call. = FALSE)

  n_blocks <- length(blocks)
  block_G <- vapply(blocks, `[[`, integer(1), "G")
  block_q <- vapply(blocks, `[[`, integer(1), "q")
  block_corr <- vapply(blocks, function(b) {
    if (.is_corr_diag(b$corr)) 0L else 1L
  }, integer(1))
  block_theta0 <- integer(n_blocks)
  theta0 <- numeric(0)
  off <- 0L
  for (b in seq_len(n_blocks)) {
    block_theta0[b] <- off
    th <- .theta0_sigma(block_q[b], blocks[[b]]$corr)
    theta0 <- c(theta0, th)
    off <- off + length(th)
  }

  # Legacy-compatible nRandComp / namesRand (one entry per SD)
  namesRand <- character(0)
  nRandComp_legacy <- integer(0)
  sigma_names <- character(0)
  for (b in seq_len(n_blocks)) {
    nb <- blocks[[b]]$name
    tn <- blocks[[b]]$term_names
    for (t in seq_len(block_q[b])) {
      sigma_names <- c(sigma_names, paste0(nb, ".", tn[t]))
    }
    # PROreg-style: one component per block with G*q effects when reporting nRandComp?
    # Keep: nRandComp = G for q=1; for q>1 store G*q in nRand for total
    nRandComp_legacy <- c(nRandComp_legacy, block_G[b] * block_q[b])
    namesRand <- c(namesRand, nb)
  }

  list(
    Z = Z,
    blocks = blocks,
    n_blocks = n_blocks,
    block_G = as.integer(block_G),
    block_q = as.integer(block_q),
    block_corr = as.integer(block_corr),
    block_theta0 = as.integer(block_theta0),
    theta0 = as.numeric(theta0),
    nRand = ncol(Z),
    namesRand = namesRand,
    sigma_names = sigma_names,
    nRandComp = as.integer(nRandComp_legacy),
    nComp = length(sigma_names),
    corr_default = corr_default
  )
}
