\name{bigmlogit}
\alias{bigmlogit}
\alias{mean.bigmlogit}
\alias{predict.bigmlogit}
\title{Multinomial logit model for big data sets}
\description{
  Estimation by maximum likelihood of the multinomial logit model, with
  alternative-specific and/or individual specific variables, for big data
  sets (millions of individuals and dozens of alternatives)
}
\usage{
bigmlogit(formula, data, chunksize = 1E3, alternatives = NULL,
         mixed = NULL, fitted = NULL, ...)
\method{mean}{bigmlogit}(x, ...)
\method{predict}{bigmlogit}(object, data, alternatives = NULL, mixed = NULL, 
    print.level = 0, ...) 
}

\arguments{
  \item{object}{a \code{bigmlogit} object,}
  \item{x}{a \code{bigmlogit} object,}
  \item{formula}{a symbolic description of the model to be estimated,}
  \item{data}{the data: a connection to a data base,}
  \item{chunksize}{the size of the blocks of data,}
  \item{alternatives}{if \code{NULL}, the alternative-specific
    covariates data frame is stored in the data base ; otherwise, it is
    provided by this argument,}
  \item{mixed}{if \code{NULL}, each mixed covariate is stored as a table
    in the data base ; otherwise, they are provided by this argument
    which is a named list containing a data frame for every mixed
    variable,}
  \item{print.level}{the amont of information printed while computing
    the predictions,}
  \item{fitted}{the name of the data base which would contain the fitted
  values,}
  \item{...}{further arguments passed to \code{maxLik}.}
}

\value{
  An object of class \code{c("bigmlogit", "maxLik")}.
}

\seealso{
  
  \code{\link{multinom}} from package \code{nnet} performs the
  estimation of the multinomial logit model with individual specific
  variables, \code{mlogit} from package \code{mlogit} with individual
  and alternative specific variables.  
}

\author{Yves Croissant}


\keyword{regression}
