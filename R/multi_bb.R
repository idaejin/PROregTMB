#' Stack outcomes for multidimensional (joint) BBmm (advanced helper)
#'
#' Low-level builder for the shared-latent multivariate **cross-sectional**
#' design. Prefer [BBmm()] with \code{cbind(y1,y2) ~ x} or long data and
#' \code{dim = "domain"}.
#'
#' @param y_list List of length \(L\); each element length-\(n\).
#' @param x Optional common covariate (if `formula_list` missing).
#' @param formula_list Optional list of \(L\) RHS formulas.
#' @param data Optional data frame.
#' @param m Scalar, length-\(L\), or length-\(n L\).
#' @param id Optional subject id (length \(n\)).
#' @return List with `X`, `y`, `Z`, `nRandComp`, `m`, `nDim`, ...
#' @export
multi_bb_stack <- function(y_list, x = NULL, formula_list = NULL,
                           data = NULL, m, id = NULL) {
  .multi_bb_stack_impl(
    y_list = y_list, x = x, formula_list = formula_list,
    data = data, m = m, id = id
  )
}

#' @keywords internal
.multi_bb_stack_impl <- function(y_list, x = NULL, formula_list = NULL,
                                 data = NULL, m, id = NULL) {
  if (!is.list(y_list) || !length(y_list)) {
    stop("y_list must be a non-empty list of response vectors", call. = FALSE)
  }
  L <- length(y_list)
  dim_names <- names(y_list)
  if (is.null(dim_names) || any(!nzchar(dim_names))) {
    dim_names <- paste0("dim", seq_len(L))
  }
  n <- length(y_list[[1]])
  if (any(vapply(y_list, length, 1L) != n)) {
    stop("All elements of y_list must have the same length n", call. = FALSE)
  }
  if (is.null(id)) {
    id <- factor(seq_len(n))
  } else {
    id <- droplevels(as.factor(id))
    if (length(id) != n) stop("id must have length n", call. = FALSE)
  }

  if (!is.null(formula_list)) {
    if (length(formula_list) != L) {
      stop("formula_list must have length nDim = length(y_list)", call. = FALSE)
    }
    if (is.null(data)) stop("data is required with formula_list", call. = FALSE)
    X_list <- vector("list", L)
    for (l in seq_len(L)) {
      f <- formula_list[[l]]
      if (length(f) == 3L) f <- f[-2L]
      mf <- stats::model.frame(f, data = data)
      Xl <- stats::model.matrix(attr(mf, "terms"), data = mf)
      if (nrow(Xl) != n) {
        stop("formula_list[[", l, "]] yields nrow != n", call. = FALSE)
      }
      colnames(Xl) <- paste0(dim_names[l], ".", colnames(Xl))
      X_list[[l]] <- Xl
    }
  } else {
    if (is.null(x)) stop("Provide x or formula_list", call. = FALSE)
    if (length(x) != n) stop("x must have length n", call. = FALSE)
    X_list <- lapply(seq_len(L), function(l) {
      Xl <- cbind("(Intercept)" = 1, x = as.numeric(x))
      colnames(Xl) <- paste0(dim_names[l], ".", colnames(Xl))
      Xl
    })
  }
  X <- as.matrix(bdiag(X_list))
  colnames(X) <- unlist(lapply(X_list, colnames))

  Z1 <- stats::model.matrix(~ id - 1)
  colnames(Z1) <- paste0("id", seq_len(ncol(Z1)))
  Z <- do.call(rbind, replicate(L, Z1, simplify = FALSE))
  y <- unlist(y_list, use.names = FALSE)
  m_obs <- .expand_m_multi(m, n = n, L = L)

  list(
    X = X, y = as.numeric(y), Z = Z, nRandComp = ncol(Z1),
    m = m_obs, nDim = as.integer(L), n = as.integer(n), id = id,
    dim_names = dim_names, beta_names = colnames(X)
  )
}

#' @keywords internal
.expand_m_multi <- function(m, n, L) {
  if (length(m) == 1L) return(rep(as.numeric(m), n * L))
  if (length(m) == L) {
    return(unlist(lapply(seq_len(L), function(l) rep(as.numeric(m)[l], n))))
  }
  if (is.list(m) && length(m) == L) {
    m_obs <- unlist(m, use.names = FALSE)
    if (length(m_obs) != n * L) stop("list m must have length n*L", call. = FALSE)
    return(as.numeric(m_obs))
  }
  if (length(m) == n * L) return(as.numeric(m))
  stop("m must be scalar, length L, length n*L, or list of L vectors",
       call. = FALSE)
}

#' @keywords internal
.is_cbind_response <- function(formula) {
  if (is.null(formula) || length(formula) < 3L) return(FALSE)
  lhs <- formula[[2L]]
  is.call(lhs) && identical(lhs[[1L]], as.name("cbind"))
}

#' @keywords internal
.resolve_m_from_data <- function(m, data, n) {
  if (is.character(m)) {
    if (length(m) != 1L || is.null(data) || !m %in% names(data)) {
      stop("m = '", paste(m, collapse = ","), "' not found in data",
           call. = FALSE)
    }
    return(as.numeric(data[[m]]))
  }
  m
}

#' @keywords internal
.rhs_formula <- function(formula) {
  tl <- attr(stats::terms(formula), "term.labels")
  if (!length(tl)) return(~ 1)
  stats::reformulate(tl)
}

#' Shared-latent multivariate design for the clean BBmm API
#' @keywords internal
.prepare_multivariate_bbmm <- function(fixed.formula, data, m, dim = NULL,
                                       random = NULL, random.formula = NULL,
                                       corr = "unstructured") {
  data <- as.data.frame(data)

  # ===== Long format: dim = "domain" =====
  if (!is.null(dim)) {
    if (!dim %in% names(data)) {
      stop("dim = '", dim, "' not found in data", call. = FALSE)
    }
    if (.is_cbind_response(fixed.formula)) {
      stop("Use either cbind(...) ~ x (wide) or dim= (long), not both",
           call. = FALSE)
    }
    dim_fac <- droplevels(as.factor(data[[dim]]))
    L <- nlevels(dim_fac)
    dim_names <- levels(dim_fac)
    ord <- order(as.integer(dim_fac))
    data_o <- data[ord, , drop = FALSE]
    dim_fac_o <- dim_fac[ord]

    y_name <- all.vars(fixed.formula[[2L]])[1]
    if (is.na(y_name) || !y_name %in% names(data_o)) {
      stop("Response in fixed.formula not found in data", call. = FALSE)
    }

    # Parametric X only; s() built later on the stacked long rows
    rhs <- .parametric_rhs(fixed.formula, drop = dim)

    m_raw <- .resolve_m_from_data(m, data_o, nrow(data_o))
    if (length(m_raw) == 1L) m_raw <- rep(as.numeric(m_raw), nrow(data_o))
    if (length(m_raw) != nrow(data_o)) {
      stop("m must be scalar, column name, or length nrow(data)", call. = FALSE)
    }

    X_list <- vector("list", L)
    y_list <- vector("list", L)
    m_list <- vector("list", L)
    for (l in seq_len(L)) {
      rows <- which(dim_fac_o == dim_names[l])
      dsub <- data_o[rows, , drop = FALSE]
      mf <- stats::model.frame(rhs, data = dsub)
      Xl <- stats::model.matrix(attr(mf, "terms"), data = mf)
      colnames(Xl) <- paste0(dim_names[l], ".", colnames(Xl))
      X_list[[l]] <- Xl
      y_list[[l]] <- as.numeric(dsub[[y_name]])
      m_list[[l]] <- as.numeric(m_raw[rows])
    }
    X <- as.matrix(bdiag(X_list))
    colnames(X) <- unlist(lapply(X_list, colnames))
    y <- unlist(y_list, use.names = FALSE)
    m_obs <- unlist(m_list, use.names = FALSE)

    if (is.null(random) && is.null(random.formula)) {
      stop("multivariate BBmm requires random = ~ (1 | id) or similar",
           call. = FALSE)
    }
    re <- .build_re_structure(
      random = random,
      random.formula = random.formula,
      data = data_o,
      corr = corr,
      nObs = length(y)
    )
    return(list(
      y = y, X = X, m = m_obs, nDim = as.integer(L),
      data_re = data_o, re = re, formula_out = fixed.formula,
      dim_names = dim_names, mode = "long"
    ))
  }

  # ===== Wide format: cbind(y1,y2) ~ x  (shared RI) =====
  if (.is_cbind_response(fixed.formula)) {
    mf <- stats::model.frame(fixed.formula, data = data)
    Y <- as.matrix(stats::model.response(mf))
    if (ncol(Y) < 2L) stop("cbind() response needs >= 2 columns", call. = FALSE)
    L <- ncol(Y)
    dim_names <- colnames(Y)
    if (is.null(dim_names) || any(!nzchar(dim_names))) {
      dim_names <- paste0("dim", seq_len(L))
      colnames(Y) <- dim_names
    }
    n <- nrow(Y)
    y_list <- lapply(seq_len(L), function(l) as.numeric(Y[, l]))
    names(y_list) <- dim_names
    # Parametric X only; s() stacked after multi prep (shared across dims)
    rhs <- .parametric_rhs(fixed.formula)
    formula_list <- replicate(L, rhs, simplify = FALSE)

    # grouping factor from random
    id <- NULL
    if (!is.null(random) && length(.findbars(random))) {
      gname <- all.vars(.findbars(random)[[1]][[3L]])
      if (length(gname) == 1L && gname %in% names(data)) id <- data[[gname]]
    }
    if (is.null(id) && !is.null(random.formula)) {
      rf <- stats::update(random.formula, ~ . - 1)
      id <- stats::model.frame(rf, data = data)[[1]]
    }
    if (is.null(id)) id <- factor(seq_len(n))

    # Wide + RI+RS: ask user to use long + dim=
    if (!is.null(random) && length(.findbars(random))) {
      lhs <- .findbars(random)[[1]][[2L]]
      if (!(identical(lhs, 1) || identical(lhs, 1L))) {
        stop(
          "Wide cbind(...) currently supports shared random intercepts ",
          "random = ~ (1 | id). For RI+RS use long data with dim = \"...\".",
          call. = FALSE
        )
      }
    }

    m_res <- .resolve_m_from_data(m, data, n)
    des <- .multi_bb_stack_impl(
      y_list = y_list, formula_list = formula_list,
      data = data, m = m_res, id = id
    )
    re <- list(
      Z = des$Z,
      nRand = ncol(des$Z),
      namesRand = "id",
      nRandComp = des$nRandComp,
      nComp = 1L,
      n_blocks = 1L,
      block_G = as.integer(des$nRandComp),
      block_q = 1L,
      block_corr = 0L,
      block_theta0 = 0L,
      theta0 = log(0.5),
      blocks = list(list(
        name = "id", G = as.integer(des$nRandComp), q = 1L, corr = "diag",
        term_names = "(Intercept)"
      )),
      sigma_names = "id.(Intercept)",
      corr_default = "diag"
    )
    return(list(
      y = des$y, X = des$X, m = des$m, nDim = des$nDim,
      data_re = data, re = re, formula_out = fixed.formula,
      dim_names = dim_names, mode = "wide"
    ))
  }

  NULL
}
