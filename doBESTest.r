BESTest <- function(A, g, N = 1e5, modelName = c('SBMpoisson', 'dcSBMpoisson'))
{
  # This is the R implementation of the method introduced by Peel et al (2017).
  # TODO: Implement all the models in the original codes

  # ----------------------------------------------------------------------------
  isValidN <- function(N)
  {
    # Check to see if N is valid
    if(N < 1)
      message('nSamples must be positive')
    else if(!(is.integer(N) | is.numeric(N)))
      message('nSamples must be integer')
    else
      T
  }

  isValidMtx <- function(A)
  {
    # Check to see if adjacency matrix is valid
    isSymmetricBool   <- isSymmetric(A)
    isUpperTriangular <- sum(A[!upper.tri(A)]) == 0
    isLowerTriangular <- sum(A[!lower.tri(A)]) == 0
    if((isSymmetricBool + isUpperTriangular + isLowerTriangular) == 0)
      message('Adjacency matrix is not symmetric or triangular')
    if(isSymmetricBool)
      return(A)
    else
      return(A+solve(A))
  }

  isValidPartition <- function(g)
  {
    # Check to see if partition is valid
    if(!(is.integer(g) | is.numeric(g)))
      message('Partition vector is non-integer. Let g(i) = integer ID of vertex i`s group.')
    k  <- min(g)
    g  <- g - k + 1
    K1 <- max(g)
    K2 <- length(unique(g))
    if(K1 != K2)
      sapply(g, function(i) which(unique(g) == g[i]))
    else
      g
  }

  SBMlikelihood <- function(A, g)
  {
    ######################
    # (dc)SBM likelihood #
    ######################
    N     <- nrow(A)
    k     <- apply(A, 1, sum)
    K     <- max(g)
    kappa <- vapply(1:K, function(i) sum(k[g == i]),   numeric(1))
    n     <- vapply(1:K, function(i) sum(  g == i) ,   numeric(1))
    theta <- vapply(1:N, function(i) k[i]/kappa[g[i]], numeric(1))
    # Omega matrix (R is clever on matrix manipulation too)
    m     <- sapply(1:K, function(i) sapply(1:K, function(j) sum(A[g==i,g==j])))

    intermediate <- m*log(diag(1/kappa)%*%m%*%diag(1/kappa))
    intermediate[is.nan(intermediate)] <- 0
    DCSBM <- sum(intermediate)

    intermediate <- m*log(diag(1/n)%*%m%*%diag(1/n))
    intermediate[is.nan(intermediate)] <- 0
    SBM <- sum(intermediate)

    return(c(SBM, DCSBM))
  }

  # ----------------------------------------------------------------------------
  model <- which(c('SBMpoisson', 'dcSBMpoisson') == modelName)
  if(!isValidN(N))
    return(NULL)
  A <- isValidMtx(A)
  g <- isValidPartition(g)

  X      <- -1*SBMlikelihood(A, g)[model]
  X.null <- sapply(1:N, function(w) -1*SBMlikelihood(A, sample(g))[model])

  return(sum(X.null >= X)/N)
}
