bigmlogit <- function(formula, data, chunksize = 1E+03, alternatives = NULL, mixed = NULL, fitted = NULL, ...){

  if (is.character(data)){
    connection <- dbConnect(dbDriver("SQLite"), dbname = data)
##    conFitted <- dbConnect(dbDriver("SQLite"), dbname = data)
  }
  else{
    connection <- data
  }
  if (!is.null(fitted)){
    conFitted <- dbConnect(dbDriver("SQLite"), dbname = fitted)
  }
  else conFitted <- NULL
  
  callT <- match.call(expand.dots = TRUE)
  formula <- Formula(formula)
  mc <- match.call(expand.dots = TRUE)
  if (is.null(mc$print.level)) plevel <- print.level <- 0 else plevel <- print.level <- mc$print.level

  # extract the model matrix for the alternative specific covariates
  if (is.null(alternatives)){
    queryZ <- dbSendQuery(connection, "select * from alternatives")
    Z <- model.matrix(formula, fetch(queryZ, n = -1), rhs = 1)[, -1, drop = FALSE]
    dbClearResult(queryZ)
  }
  else Z <- model.matrix(formula, alternatives, rhs = 1)[, -1, drop = FALSE]
  KZ <- ncol(Z)
  J <- nrow(Z)

  # for the "mixed variables", extract the names of the variables from
  # the formula and then fetch all the observations from the
  # corresponding tables and store these variables in a list
  MNames <- attr(terms(formula(formula, rhs = 3, lhs = 0)), "term.labels")
  KM <- length(MNames)
  M <- vector(length = KM, mode = "list")
  for (k in 1:KM){
    if (is.null(mixed)){
      queryM <- dbSendQuery(connection, paste("select * from", MNames[k]))
      M[[k]] <- fetch(queryM, n = - 1)
      dbClearResult(queryM)
    }
    else M[[k]] <- mixed[[MNames[[k]]]]
  }
  names(M) <- MNames
  
  # to get starting values and the names of the coefficients, estimate
  # the model on a single chunk
  queryX <- dbSendQuery(connection, "select * from individuals")
  dataX <- fetch(queryX, n = chunksize)
  responseName <- as.character(attr(formula, "lhs")[[1]]) 
  y <- dataX[[responseName]]
  X <- model.matrix(formula, dataX, rhs = 2)
  KX <- ncol(X)
  dbClearResult(queryX)

  # names of the covariates
  namesCovariates <- c(paste(rep(colnames(X), KZ), rep(colnames(Z), each = KX), sep = ":"), names(M))

  # cbind the X matrix KZ times
  X <- matrix(rep(X, KZ), nrow(X))
  # create a list of J matrices of covariates : the covariates
  # consist on every individual specific variable times every
  # alternative specific variables
  X <- lapply(seq_len(nrow(Z)),
              function(i) rep(Z[i, ], each = chunksize * KX) * X)
  # create a list containing lists of J vectors of "mixed variables"
  Mmat <- model.matrix(formula, dataX, rhs = 3)[, - 1, drop = FALSE]
  MI <- vector(length = KM, mode = "list")
  for (k in 1:KM) MI[[k]] <- lapply(M[[k]], function(x) x[Mmat[, k, drop = FALSE]])
  # cbind these variables to the other covariates
  for (k in 1:KM) X <- lapply(seq_len(J),
                              function(i) cbind(X[[i]], MI[[k]][[i]]))

  lnlstart <- function(param, X, y){
    eXb <- lapply(X, function(x) as.numeric(exp(crossprod(t(x), param))))
    seXb <- Reduce("+", eXb)
    P <- sapply(eXb, function(x) x / seXb)
    Y <- sapply(seq_len(J), function(i) y == i)
    Pch <- Reduce("+", lapply(seq_len(J), function(i) Y[, i] * P[, i]))
    lnl <-  sum(log(Pch))
    PX <- Reduce("+", lapply(seq_len(J), function(i) X[[i]] * P[, i]))
    Xch <- Reduce("+", lapply(seq_len(J), function(i) X[[i]] * Y[, i]))
    gradient <-  apply(Xch - PX, 2, sum)
    XmPX <- lapply(X, function(x) x - PX)
    hessian <- - Reduce("+", lapply(seq_len(J),
                                    function(i) crossprod(P[, i] * XmPX[[i]], XmPX[[i]])))
    attr(lnl, "gradient") <-  gradient
    attr(lnl, "hessian") <-   hessian
    lnl
  }

  starting.values <- rep(0, ncol(X[[1]]))
  names(starting.values) <- namesCovariates

  if (plevel > 1) cat("Computing the starting values estimating the model on one chunk\n")
  result <- maxLik(lnlstart, start = starting.values, X = X, y = y, ...)
  result <- maxLik(lnlSQLite, start = coef(result), chunksize = chunksize,
                   connection = connection, conFitted = conFitted, Z = Z, M = M,
                   formula = formula, plevel = print.level, ...)
  dbDisconnect(connection)
  if (!is.null(fitted)){
    dbDisconnect(conFitted)
  }
  result$call <- callT
  result$alternatives <- Z
  result$mixed <- M
  result$chunksize <- chunksize
  result$formula <- formula
  result$fitted <- attr(result$maximum, "fitted")
  result$mean <- attr(result$maximum, "mean")
  class(result) <- c("bigmlogit", "maxLik", "maxim")
  result
}

fitted.bigmlogit <- function(object, ...){
  object$fitted
}

predict.bigmlogit <- function(object, data, alternatives = NULL, mixed = NULL, print.level = 0, ...){
  if (is.character(data)){
    data <- dbConnect(dbDriver("SQLite"), dbname = data)
  }
  formula <- object$formula
  chunksize <- object$chunksize
  if (is.null(alternatives) & is.null(mixed)) object$fitted
  else{
    if (is.null(alternatives)) Z <- object$alternatives
    else Z <- model.matrix(formula, alternatives, rhs = 1)[, -1, drop = FALSE]
    if (is.null(mixed)) M <- object$mixed
    else M <- mixed
    z <- object$chunksize
    result <- lnlSQLite(coef(object), chunksize = chunksize, connection = data,
                        conFitted = NULL, 
                        Z = Z, M = M, formula = formula, plevel = print.level)
    result <- attr(result, "fitted")
  }
  dbDisconnect(data)
  result
}


mean.bigmlogit <- function(x, ...) x$mean
