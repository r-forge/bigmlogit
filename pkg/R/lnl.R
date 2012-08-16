lnlSQLite <- function(param, chunksize = 1E+03, connection, conFitted, Z, M, formula, plevel){
  # param is the vector of parameters, size the size of the chunks,
  # connection the connection to the data base and formula a Formula
  # object which describes the model to be estimated
  
  # get the number of individuals and compute the number of chunks
  querySize <- dbSendQuery(connection, "select count(*) from individuals")
  n <- as.numeric(fetch(querySize))
  dbClearResult(querySize)
  L <- ceiling(n / chunksize)
  
  # for the individual-specific variables, just write the query ; the
  # data will later be fetched by chunks of size defined by the size
  # argument
  queryX <- dbSendQuery(connection, "select * from individuals")

  # the lnl, gradient and hessian are updated at each iteration simply
  # by adding the contribution of the new chunk ; for now, just define
  # them as 0
  gradient <- hessian <- lnl <- 0
  fvalues <- ymean <- 0
  
  stillToLoad <- n
  if (plevel > 1) cat('loading data : ')
  for (l in 1:L){
    # the size of the chunk is chunksize except for the last one
    if(stillToLoad < chunksize) chunksize <- stillToLoad
    prct <- round( (n - stillToLoad) / n * 100, 2)
    if (plevel > 1) cat(prct, '% ')
    # get the data frame for a chunk
    dataX <- fetch(queryX, n = chunksize)
    # create the model matrix
    Xmatrix <- model.matrix(formula, dataX, rhs = 2)
    namesX <- colnames(Xmatrix)
    KX <- ncol(Xmatrix)
    # extract the response from the data frame
    responseName <- as.character(attr(formula, "lhs")[[1]]) 
    y <- dataX[[responseName]]
    # cbind the individual-specific variables matrix KZ times
    X <- c()
    KZ <- ncol(Z)
    J <- nrow(Z)
    for (z in 1:KZ) X <- cbind(X, Xmatrix)
    # create a list of J matrices of covariates : the covariates
    # consist on every individual specific variable times every
    # alternative specific variables
    X <- lapply(seq_len(nrow(Z)),
                function(i) rep(Z[i, ], each = chunksize * KX) * X)
    # create a list containing lists of J vectors of "mixed variables"
    Mmat <- model.matrix(formula, dataX, rhs = 3)[, - 1, drop = FALSE]
    KM <- length(M)
    mixedVariables <- vector(length = KM, mode = "list")
    for (k in 1:KM) mixedVariables[[k]] <- lapply(M[[k]],
                                                  function(x) x[Mmat[, k, drop = FALSE]])
    # cbind these variables to the other covariates
    for (k in 1:KM) X <- lapply(seq_len(J),
                                function(i) cbind(X[[i]], mixedVariables[[k]][[i]]))
    # create a list of length J containing the exp of the level of the
    # deterministic part of utility for every individual (rows)
    eXb <- lapply(X, function(x) as.numeric(exp(crossprod(t(x), param))))
    # compute the denominator of the logit probabilities
    seXb <- Reduce("+", eXb)
    # compute then a list of length J containing the probabilities
    P <- sapply(eXb, function(x) x / seXb)
    # append the probabilities for this chunk on the Fitted table
    if (!is.null(conFitted)){
      dbWriteTable(conFitted, "fitted", data.frame(P), append = TRUE, row.names=FALSE)
    }
    Pd <- as.data.frame(P)
    fvalues <- fvalues + apply(P, 2, sum)
    # create a list of length J containing logical vectors indicating
    # the chosen alternative
    Y <- sapply(seq_len(J), function(i) y == i)
    ymean <- ymean + apply(Y, 2, sum)
    # compute a vector containing the probability for the chosen
    # alternative
    Pch <- Reduce("+", lapply(seq_len(J), function(i) Y[, i] * P[, i]))
    # the log likelihood is simply the sum of the log of these
    # probabilities ; update it with the contribution of the present
    # chunk
    lnl <- lnl + sum(log(Pch))
    # to compute the gradient, on needs the matrix of the covariates
    # averaged by the probabilities of every alternative
    PX <- Reduce("+", lapply(seq_len(J), function(i) X[[i]] * P[, i]))
    # to compute the gradient, on needs the matrix of the covariates
    # for the chosen alternative
    Xch <- Reduce("+", lapply(seq_len(J), function(i) X[[i]] * Y[, i]))
    # update the gradient with the contribution of the present chunk
    gradient <-  gradient + apply(Xch - PX, 2, sum)
    # update the hessian with the contribution of the present chunk
    XmPX <- lapply(X, function(x) x - PX)
    hessian <- hessian -
      Reduce("+", lapply(seq_len(J),
                         function(i) crossprod(P[, i] * XmPX[[i]], XmPX[[i]])))
    stillToLoad <- stillToLoad - chunksize
  }
  if (plevel > 1) cat('\n')
  # disconnect before exit
  dbClearResult(queryX)
  # return the log-likelihood, with the gradient and the hessian as
  # attributes
  attr(lnl, "gradient") <-  gradient
  attr(lnl, "hessian") <-  hessian
  attr(lnl, "fitted") <- fvalues / n
  attr(lnl, "mean") <- ymean / n
  lnl
}

