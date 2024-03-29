#' Collective Matrix Factorization (CMF)
#'
#' Collective matrix factorization (CMF) finds joint low-rank
#' representations for a collection of matrices with shared
#' row or column entities. This package learns a variational
#' Bayesian approximation for CMF, supporting multiple
#' likelihood potentials and missing data, while identifying
#' both factors shared by multiple matrices and factors
#' private for each matrix.
#'
#' This package implements a variational Bayesian approximation for
#' CMF, following the presentation in "Group-sparse embeddings in
#' collective matrix factorization" (see references below).
#'
#' The main functionality is provided by the function
#' [CMF()] that is used for learning the model, and by the
#' function [predictCMF()] that estimates missing entries
#' based on the learned model. These functions take as input
#' lists of matrices in a specific sparse format that stores
#' only the observed entries but that explicitly stores
#' zeroes (unlike most sparse matrix representations).
#' For converting between regular matrices and this sparse
#' format see [matrix_to_triplets()] and [triplets_to_matrix()].
#'
#' The package can also be used to learn Bayesian canonical
#' correlation analysis (CCA) and group factor analysis (GFA)
#' models, both of which are special cases of CMF. This is likely to be
#' useful for people looking for CCA and GFA solutions supporting
#' missing data and non-Gaussian likelihoods.
#'
#' @author Arto Klami \email{arto.klami@@cs.helsinki.fi} and Lauri Väre
#'
#' Maintainer: Felix Held \email{felix.held@@gmail.se}
#'
#' @references
#' Arto Klami, Guillaume Bouchard, and Abhishek Tripathi.
#' Group-sparse embeddings in collective matrix factorization.
#' arXiv:1312.5921, 2013.
#'
#' Arto Klami, Seppo Virtanen, and Samuel Kaski.
#' Bayesian canonical correlation analysis. Journal of Machine
#' Learning Research, 14(1):965--1003, 2013.
#'
#' Seppo Virtanen, Arto Klami, Suleiman A. Khan, and Samuel Kaski.
#' Bayesian group factor analysis. In Proceedings of the 15th
#' International Conference on Artificial Intelligence and Statistics,
#' volume 22 of JMLR:W&CP, pages 1269-1277, 2012.
#'
#' @examples
#' require("CMF")
#'
#' # Create data for a circular setup with three matrices and three
#' # object sets of varying sizes.
#' X <- list()
#' D <- c(10, 20, 30)
#' inds <- matrix(0, nrow = 3, ncol = 2)
#'
#' # Matrix 1 is between sets 1 and 2 and has continuous data
#' inds[1, ] <- c(1, 2)
#' X[[1]] <- matrix(
#'   rnorm(D[inds[1, 1]] * D[inds[1, 2]], 0, 1),
#'   nrow = D[inds[1, 1]]
#' )
#'
#' # Matrix 2 is between sets 1 and 3 and has binary data
#' inds[2, ] <- c(1, 3)
#' X[[2]] <- matrix(
#'   round(runif(D[inds[2, 1]] * D[inds[2, 2]], 0, 1)),
#'   nrow = D[inds[2, 1]]
#' )
#'
#' # Matrix 3 is between sets 2 and 3 and has count data
#' inds[3, ] <- c(2, 3)
#' X[[3]] <- matrix(
#'   round(runif(D[inds[3, 1]] * D[inds[3, 2]], 0, 6)),
#'   nrow = D[inds[3, 1]]
#' )
#'
#' # Convert the data into the right format
#' triplets <- lapply(X, matrix_to_triplets)
#'
#' # Missing entries correspond to missing rows in the triple representation
#' # so they can be removed from training data by simply taking a subset
#' # of the rows.
#' train <- list()
#' test <- list()
#' keep_for_training <- c(100, 200, 300)
#' for (m in 1:3) {
#'   subset <- sample(nrow(triplets[[m]]), keep_for_training[m])
#'   train[[m]] <- triplets[[m]][subset, ]
#'   test[[m]] <- triplets[[m]][-subset, ]
#' }
#'
#' # Learn the model with the correct likelihoods
#' K <- 4
#' likelihood <- c("gaussian", "bernoulli", "poisson")
#' opts <- getCMFopts()
#' opts$iter.max <- 500 # Less iterations for faster computation
#' model <- CMF(train, inds, K, likelihood, D, test = test, opts = opts)
#'
#' # Check the predictions
#' # Note that the data created here has no low-rank structure,
#' # so we should not expect good accuracy.
#' print(test[[1]][1:3, ])
#' print(model$out[[1]][1:3, ])
#'
#' # predictions for the test set using the previously learned model
#' out <- predictCMF(test, model)
#' print(out$out[[1]][1:3, ])
#' print(out$error[[1]])
#' # ...this should be the same as the output provided by CMF()
#' print(model$out[[1]][1:3, ])
#'
#' @docType package
#' @name CMF-package
#' @useDynLib CMF, .registration = TRUE
NULL

#' Default options for CMF
#'
#' A helper function that creates a list of options to be passed to `CMF`.
#' To run the code with other option values, first run this function and
#' then directly modify the entries before passing the list to `CMF`.
#'
#' Most of the parameters are for controlling the optimization, but some will
#' alter the model itself. In particular, `useBias` is used for turning
#' the bias terms on and off, and `method` will change the prior for `U`.
#'
#' The default choice for `method` is `"gCMF"`, providing the
#' group-wise sparse CMF that identifies both shared and private factors
#' (see Klami et al. (2013) for details). The value `"CMF"` turns off
#' the group-wise sparsity, providing a CMF solution that attempts to learn
#' only factors shared by all matrices. Finally, `method="GFA"` implements
#' the group factor analysis (GFA) method, by fixing the variance of
#' `U[[1]]` to one and forcing `useBias=FALSE`. Then `U[[1]]` can be
#' interpreted as latent variables with unit variance and zero mean,
#' as assumed by GFA and CCA (special case of GFA with `M = 2`). Note that as a
#' multi-view learning method `"GFA"` requires all matrices to share the
#' same rows, the very first entity set.
#'
#' @return Returns a list of:
#' \item{init.tau}{Initial value for the noise precisions. Only matters for
#'                 Gaussian likelihood.}
#' \item{init.alpha}{Initial value for the automatic relevance determination
#'                   (ARD) prior precisions.}
#' \item{grad.reg }{The regularization parameter for the under-relaxed Newton
#'                  iterations. 0 = no regularization, larger values provide
#'                  increasing regularization. The value must be below 1.}
#' \item{gradIter}{How many gradient steps for updating the projections are
#'                  performed during each iteration of the whole algorithm.
#'                  Default is 1.}
#' \item{grad.max}{Maximum absolute change for the elements of the projection
#'                 matrices during one gradient step. Small values help to
#'                 prevent over-shooting, wheres inf results to no constraints.
#'                 Default is `inf`.}
#' \item{iter.max }{Number of iterations for the whole algorithm.}
#' \item{computeCost}{Should the cost function values be computed or not.
#'                    Defaults to `TRUE`.}
#' \item{verbose}{0 = supress all printing, 1 = print current iteration and
#'                test RMSE every now and then, 2 = in addition to level 1
#'                print also the current gradient norm.}
#' \item{useBias}{Set this to `FALSE` to exclude the row and column bias terms.
#'                The default is `TRUE`.}
#' \item{method}{Default value of "gCMF" computes the CMF with group-sparsity.
#'               The other possible values are "CMF" for turning off the
#'               group-sparsity prior, and "GFA" for implementing group factor
#'               analysis (and canonical correlation analysis when `M = 2`).}
#' \item{prior.alpha_0}{Hyperprior values for the gamma prior for ARD.}
#' \item{prior.alpha_0t}{Hyperprior values for the gamma prior for tau.}
#' @author Arto Klami and Lauri Väre
#' @seealso 'CMF'
#' @references
#' Arto Klami, Guillaume Bouchard, and Abhishek Tripathi.
#' Group-sparse embeddings in collective matrix factorization.
#' arXiv:1312.5921, 2014.
#'
#' Seppo Virtanen, Arto Klami, Suleiman A. Khan, and Samuel Kaski.
#' Bayesian group factor analysis. In Proceedings of the 15th
#' International Conference on Artificial Intelligence and Statistics,
#' volume 22 of JMLR:W&CP, pages 1269-1277, 2012.
#' @examples
#'
#' CMF_options <- getCMFopts()
#' CMF_options$iter.max <- 500 # Change the number of iterations from default
#'                             # of 200 to 500.
#' CMF_options$useBias <- FALSE # Do not take row and column means into
#'                              # consideration.
#' # These options will be in effect when CMF_options is passed on to CMF.
#'
#' @export
getCMFopts <- function() {
  # Initial value for the noise precisions (only matters for
  # Gaussian likelihood)
  init.tau <- 10

  # Initial value for the prior precisions
  init.alpha <- 5

  # Parameters for the gradient optimization of the projection
  # matrices.
  # The algorithm is Newton-Raphson with diagonal Hessian
  # and successive under-relaxation (= Richardson extrapolation),
  # so that
  #   theta_new = grad.reg*theta_old + (1-grad.reg)(theta_old - g/h)
  # where g is gradient and h is Hessian. grad.reg should be between 0 and 1
  # to stabilize the algorithm. 0 overshoots, 1 does not update anything.
  #
  # If grad.max < Inf, the step lengths are capped at grad.Max
  # grad.iter tells how many gradient steps to take within each update.
  grad.reg <- 0.7
  grad.iter <- 1
  grad.max <- Inf

  # Parameters for controlling when the algorithm stops.
  iter.max <- 200

  computeCost <- TRUE
  verbose <- 1 # 1=print progress every now and then, 0=supress all printing
  useBias <- TRUE # Whether to include bias terms;
                  # useBias = FALSE means no bias terms

  # Whether to do GFA/CCA instead of CMF
  # If this is set to "GFA" then useBias=FALSE and alpha is set to one for
  # the first entity set, so that U[[1]] become latent variables with zero mean
  # and unit variance.
  # It this is set fo "CMF" the alpha-parameter becomes the same for all
  # matrices.
  method <- "gCMF"

  # Hyperparameters
  # - alpha_0, beta_0 for the ARD precisions
  # - alpha_0t, beta_0t for the residual noise predictions
  prior.alpha_0 <- prior.beta_0 <- 1
  prior.alpha_0t <- prior.beta_0t <- 0.001

  return(list(
    init.tau = init.tau,
    init.alpha = init.alpha,
    grad.iter = grad.iter,
    iter.max = iter.max,
    verbose = verbose,
    computeCost = computeCost,
    grad.reg = grad.reg,
    grad.max = grad.max,
    useBias = useBias,
    method = method,
    prior.alpha_0 = prior.alpha_0,
    prior.beta_0 = prior.beta_0,
    prior.alpha_0t = prior.alpha_0t,
    prior.beta_0t = prior.beta_0t
  ))
}


#' Predict with CMF
#'
#' Code for predicting missing elements with an existing CMF model.
#' The predictions are made for all of the elements specified in the list of
#' input matrices `X`. The function also returns the root mean square error
#' (RMSE) between the predicted outputs and the values provided in `X`.
#'
#' Note that `X` needs to be provided as a set of triplets instead of as
#' a regular matrix. See [matrix_to_triplets()].
#'
#' @param X A list of sparse matrices specifying the indices for which to
#'          make the predictions.
#'          These matrices must correspond to the structure used for `X`
#'          when learning the model with `CMF`.
#' @param model A list of model parameter values provided by `CMF`.
#' @return A list of
#' \item{out}{A list of matrices corresponding to predictions for each
#'            matrix in `X`.}
#' \item{error}{A vector containing the root-mean-square error for each
#'              matrix separately.}
#'
#' @author Arto Klami and Lauri Väre
#' @examples
#'
#' # See CMF-package for an example.
#'
#' @export
predictCMF <- function(X, model) {
  D <- model$D
  inds <- model$inds

  for (i in seq_along(X)) {
    if(!p_check_sparsity(
      X[[i]], D[inds[i, 1]], D[inds[i, 2]])) {
        stop(paste0("Input matrix ", i, " is not in the correct format."))
    }
  }

  out <- list()
  for (m in 1:model$M) {
    indices <- matrix(as.integer(X[[m]][, 1:2]), ncol = 2)
    xi <- p_updatePseudoData(
      indices, model$U[[inds[m, 1]]], model$U[[inds[m, 2]]],
      model$bias[[m]]$row$mu, model$bias[[m]]$col$mu
    )
    out[[m]] <- unname(cbind(indices, xi))
  }
  for (m in which(model$likelihood == "bernoulli")) {
    out[[m]][, 3] <- exp(out[[m]][, 3])
    out[[m]][, 3] <- out[[m]][, 3] / (1 + out[[m]][, 3])
  }
  for (m in which(model$likelihood == "poisson")) {
    out[[m]][, 3] <- exp(out[[m]][, 3])
    out[[m]][, 3] <- log(1 + out[[m]][, 3])
  }

  error <- rep(0, length(X))
  for (m in seq_along(X)) {
    for (r in seq_len(nrow(X[[m]]))) {
      error[m] <- error[m] + (X[[m]][r, 3] - out[[m]][r, 3])^2
    }
    error[m] <- sqrt(error[m] / nrow(X[[m]]))
  }

  return(list(out = out, error = error))
}

#' Collective Matrix Factorization
#'
#' Learns the CMF model for a given collection of M matrices.
#' The code learns the parameters of a variational approximation for CMF,
#' and also computes predictions for indices specified in `test`.
#'
#' The variational approximation is fully factorized over all of the model
#' parameters, including individual elements of the projection matrices.
#' The parameters for the projection matrices are updated jointly by
#' Newton-Raphson method, whereas the rest use closed-form updates.
#'
#' Note that the input data needs to be given in a specific sparse format.
#' See [matrix_to_triplets()] for details.
#'
#' The behavior of the algorithm can be modified via the `opts` parameter.
#' See [getCMFopts()] for details. Of particular interest are the elements
#' `useBias` and `method`.
#'
#' For full description of the output parameters, see the referred publication.
#' The notation in the code follows roughly the notation used in the paper.
#'
#' @param X List of input matrices.
#' @param inds A `length(X)` times 2 matrix that links dimensions of the
#'             matrices in `X` to object sets. `inds[m, 1]` tells which
#'             object set corresponds to the rows in matrix `X[[m]]`,
#'             and `inds[m, 2]` tells the same for the columns.
#' @param K The number of factors.
#' @param opts A list of options as given by [getCMFopts()].
#'             If set to `NULL`, the default values will be used.
#' @param likelihood A list of likelihood choices, one for each matrix in X.
#'                   Each entry should be a string with possible values of:
#'                   "gaussian", "bernoulli" or "poisson".
#' @param D A vector containing sizes of each object set.
#' @param test A list of test matrices. If not NULL, the code will compute
#'             predictions for these elements of the matrices. This duplicates
#'             the functionality of [predictCMF()].
#' @return A list of
#' \item{U}{A list of the mean parameters for the rank-K projection matrices,
#'           one for each object set.}
#' \item{covU}{A list of the variance parameters for the rank-K projection
#'             matrices, one for each object set.}
#' \item{tau}{A vector of the precision parameter means.}
#' \item{alpha}{A vector of the ARD parameter means.}
#' \item{cost}{A vector of variational lower bound values.}
#' \item{inds}{The input parameter `inds` stored for further use.}
#' \item{errors}{A vector containing root-mean-square errors for each
#'               iteration, computed over the elements indicated by the
#'               `test` parameter.}
#' \item{bias}{A list (of lists) storing the parameters of the row and
#'             column bias terms.}
#' \item{D}{The sizes of the object sets as given in the parameters.}
#' \item{K}{The number of components as given in the parameters.}
#' \item{Uall}{Matrices of U joined into one sum(D) by K matrix, for
#'             easier plotting of the results.}
#' \item{items}{A list containing the running number for each item among
#'              all object sets. This corresponds to rows of the `Uall`
#'              matrix. Each part of the list contains a vector that has the
#'              numbers for each particular object set.}
#' \item{out}{If test matrices were provided, returns the reconstructed data
#'            sets. Otherwise returns `NULL`.}
#' \item{M}{The number of input matrices.}
#' \item{likelihood}{The likelihoods of the matrices.}
#' \item{opts}{The options used for running the code.}
#' @author Arto Klami and Lauri Väre
#' @references
#' Arto Klami, Guillaume Bouchard, and Abhishek Tripathi.
#' Group-sparse embeddings in collective matrix factorization.
#' arXiv:1312.5921, 2014.
#' @examples
#' # See CMF-package for an example.
#'
#' @importFrom stats rnorm
#' @export
CMF <- function(X, inds, K, likelihood, D, test = NULL, opts = NULL) {
  if (is.null(opts)) {
    opts <- getCMFopts()
  }

  for (i in seq_along(X)) {
    if (!p_check_sparsity(X[[i]], D[inds[i, 1]], D[inds[i, 2]])) {
      stop(paste0("Input matrix ", i, " is not in the correct format."))
    }
  }

  # Store indices separately as integers
  indices <- lapply(X, function(x) matrix(as.integer(x[, 1:2]), ncol = 2))

  if (opts$method == "GFA") {
    opts$useBias <- FALSE
    if (!all(inds[, 1] == 1)) {
      stop(paste0(
        "GFA requires all matrices to share the first entity set, ",
        "since it is a multi-view learning method."
      ))
    }
  }

  # Construct the index sets and the symmetric matrix
  #
  # X are lists of observed values
  # M is the number of data sets
  # C is the number of itemsets
  # D[i] is the size of the c-th itemset
  # items[[inds[m, 1]]] tells which samples correspond to rows in X[m]
  # and items[[inds[m, 2]]] is the same for columns
  M <- length(X)
  indvec <- as.vector(inds)
  types <- seq_len(max(indvec))
  C <- max(types)
  items <- list()

  N <- 0
  for (t in types) {
    count <- D[t]
    items[[t]] <- N + seq_len(count)
    N <- N + count
  }

  alpha_0 <- opts$prior.alpha_0   # Easier access for hyperprior values
  beta_0 <- opts$prior.beta_0
  alpha_0t <- opts$prior.alpha_0t
  beta_0t <- opts$prior.beta_0t

  # ARD and noise parameters
  alpha <- matrix(opts$init.alpha, C, K)  # The mean of the ARD precisions
  b_ard <- matrix(0, C, K)                # The parameters of the Gamma
                                          #   distribution
  a_ard <- alpha_0 + D / 2                # for ARD precisions
  tau <- rep(opts$init.tau, M)            # The mean noise precisions
  a_tau <- rep(alpha_0t, M)               # The parameters of the Gamma
                                          #   distribution
  for (m in seq_len(M)) {
     a_tau[m] <- a_tau[m] + nrow(X[[m]]) / 2
  }
  b_tau <- rep(0, M)                      #     for the noise precisions

  # The projections
  U <- vector("list", length = C)
  covU <- vector("list", length = C)      # The covariances

  bias <- vector("list", length = M)
  for (i in 1:C) {
    # Random initialization
    U[[i]] <- matrix(rnorm(D[i] * K, 0, 1 / sqrt(opts$init.alpha)), D[i], K)
    covU[[i]] <- matrix(1 / D[i], D[i], K)
  }

  #
  # Parameters for modeling the row and column bias terms
  #
  for (m in 1:M) {
    bias[[m]]$row$mu <- rep(0, D[inds[m, 1]])
    bias[[m]]$row$m <- 0
    if (!opts$useBias) {
      bias[[m]]$row$nu <- rep(0, D[inds[m, 1]])
      bias[[m]]$row$lambda <- 0
      bias[[m]]$row$scale <- 0
    } else {
      bias[[m]]$row$nu <- rep(1, D[inds[m, 1]])
      bias[[m]]$row$lambda <- 1
      bias[[m]]$row$scale <- 1
    }

    bias[[m]]$col$mu <- rep(0, D[inds[m, 2]])
    bias[[m]]$col$m <- 0
    if (!opts$useBias) {
      bias[[m]]$col$nu <- rep(0, D[inds[m, 2]])
      bias[[m]]$col$lambda <- 0
      bias[[m]]$col$scale <- 0
    } else {
      bias[[m]]$col$nu <- rep(1, D[inds[m, 2]])
      bias[[m]]$col$lambda <- 1
      bias[[m]]$col$scale <- 1
    }
  }

  # Pseudo data
  origX <- lapply(X, function(x) x[, 3])
  for (m in which(likelihood == "bernoulli")) {
    tau[m] <- 0.25
    xi <- p_updatePseudoData(
      indices[[m]], U[[inds[m, 1]]], U[[inds[m, 2]]],
      bias[[m]]$row$mu, bias[[m]]$col$mu
    )
    X[[m]][, 3] <- (
      xi - (1 / (1 + exp(-xi)) - origX[[m]]) / tau[m]
    )
  }

  for (m in which(likelihood == "poisson")) {
    maxval <- max(X[[m]][, 3])
    # print(paste("Maximal count in view ",m," is ",maxval,
    #             "; consider clipping if this is very large.",sep=""))
    tau[m] <- 0.25 + 0.17 * maxval
    xi <- p_updatePseudoData(
      indices[[m]], U[[inds[m, 1]]], U[[inds[m, 2]]],
      bias[[m]]$row$mu, bias[[m]]$col$mu
    )
    X[[m]][, 3] <- (
      xi - 1 / (1 + exp(-xi)) * (1 - origX[[m]] / log(1 + exp(xi))) / tau[m]
    )
  }

  cost <- vector()  # For storing the lower bounds

  #
  # The main loop
  #
  errors <- array(0, c(opts$iter.max, M))
  for (iter in seq_len(opts$iter.max)) {
    #gc() # Just to make sure there are no memory leaks
    if ((opts$verbose > 0) && (iter %% 10 == 1)) {
      cat(paste0("Iteration: ", iter, "/", opts$iter.max, "\n"))
    }
    #
    # Update the weight matrices, using gradients (and diagonal Hessian)
    #
    norm <- 0
    for (i in sample(C)) {
      #
      # Update the variance
      #
      covU[[i]] <- matrix(alpha[i, ], nrow = D[i], ncol = K, byrow = TRUE)
      for (m in which(inds[, 1] == i)) {
        # NOTE: Directly modifies covU[[i]]
        p_covUsparse(
          X[[m]], covU[[i]], U[[inds[m, 2]]], covU[[inds[m, 2]]], 1, tau[m]
        )
      }
      for (m in which(inds[, 2] == i)) {
        p_covUsparse(
          X[[m]], covU[[i]], U[[inds[m, 1]]], covU[[inds[m, 1]]], 2, tau[m]
        )
      }
      covU[[i]] <- 1 / covU[[i]]
      # Update U itself
      par <- list(
        D = D, alpha = alpha, tau = tau,
        X = X, U = U, covU = covU, inds = inds, bias = bias
      )
      par$this <- i

      # The inverse Hessian happens to be the covariance
      # (see Ilin & Raiko, JMLR 2010)
      scale <- covU[[i]]
      reg <- opts$grad.reg
      for (n in 1:opts$grad.iter) {
        g <- p_gradUsparseWrapper(U[[i]], par, stochastic = FALSE)
        if (max(abs(scale * g)) > opts$grad.max) {
          scale <- opts$grad.max * scale / max(abs(scale * g))
        }

        U[[i]] <- reg * U[[i]] + (1 - reg) * (U[[i]] - scale * g)
        norm <- norm + sum(g^2)
      }
    }
    if ((opts$verbose > 1) && (iter %% 10 == 1)) {
      cat(paste0(" Gradient norm: ", norm, "\n"))
    }

    #
    # Update the mean profiles (and their priors)
    #
    if (opts$useBias) {
      for (m in seq_len(M)) {
        # The approximations for each data point
        temp <- p_updateMean(
          X[[m]], U[[inds[m, 1]]], U[[inds[m, 2]]], 1, bias[[m]]$col$mu
        )
        bias[[m]]$row$count <- temp$count
        bias[[m]]$row$nu <- (
          1 / (tau[m] * bias[[m]]$row$count + 1 / bias[[m]]$row$scale)
        )
        bias[[m]]$row$mu <- (
          (tau[m] * temp$sum + bias[[m]]$row$m / bias[[m]]$row$scale)
          * bias[[m]]$row$nu
        )
        # The approximation for the mean
        bias[[m]]$row$lambda <- (
          1 / (1 + 1 / bias[[m]]$row$scale * D[inds[m, 1]])
        )
        bias[[m]]$row$m <- (
          1 / bias[[m]]$row$scale * sum(bias[[m]]$row$mu) * bias[[m]]$row$lambda
        )

        # The scale
        bias[[m]]$row$scale <- (
          mean(bias[[m]]$row$nu + (bias[[m]]$row$mu - bias[[m]]$row$m)^2)
        )

        # The approximations for each data point
        temp <- p_updateMean(
          X[[m]], U[[inds[m, 1]]], U[[inds[m, 2]]], 2, bias[[m]]$row$mu
        )
        bias[[m]]$col$count <- temp$count
        bias[[m]]$col$nu <- (
          1 / (tau[m] * bias[[m]]$col$count + 1 / bias[[m]]$col$scale)
        )
        bias[[m]]$col$mu <- (
          (tau[m] * temp$sum + bias[[m]]$col$m / bias[[m]]$col$scale)
          * bias[[m]]$col$nu
        )

        # The approximation for the mean
        bias[[m]]$col$lambda <- (
          1 / (1 + 1 / bias[[m]]$col$scale * D[inds[m, 2]])
        )
        bias[[m]]$col$m <- (
          1 / bias[[m]]$col$scale * sum(bias[[m]]$col$mu) * bias[[m]]$col$lambda
        )

        # The scale
        bias[[m]]$col$scale <- (
          mean(bias[[m]]$col$nu + (bias[[m]]$col$mu - bias[[m]]$col$m)^2)
        )
      }
    }

    #
    # Update alpha, the ARD parameters
    #
    for (i in seq_len(C)) {
      tmp <- rep(beta_0 * 2, K)
      for (k in seq_len(K)) {
        tmp[k] <- tmp[k] + sum(U[[i]][, k]^2) + sum(covU[[i]][, k])
      }
      tmp <- tmp / 2
      alpha[i, ] <- a_ard[i] / tmp
      b_ard[i, ] <- tmp
    }
    if (opts$method == "CMF") {
      for (k in seq_len(K)) {
        alpha[, k] <- sum(a_ard) / sum(b_ard[, k])
      }
    }
    if (opts$method == "GFA") {
      alpha[1, ] <- 1
    }

    #
    # Update tau, the noise precisions; only needed for Gaussian likelihood
    #
    for (m in which(likelihood == "gaussian")) {
      b_tau[m] <- beta_0t
      v1 <- inds[m, 1]
      v2 <- inds[m, 2]
      b_tau[m] <- beta_0t + 0.5 * p_updateTau(
        X[[m]], U[[v1]], U[[v2]],
        covU[[v1]], covU[[v2]],
        bias[[m]]$row$mu, bias[[m]]$col$mu,
        bias[[m]]$row$nu, bias[[m]]$col$nu
      )
      tau[m] <- a_tau[m] / b_tau[m]
    }

    #
    # Update the pseudo-data; only needed for non-Gaussian likelihoods
    #
    #print("Update pseudo data")
    for (m in which(likelihood == "bernoulli")) {
      xi <- p_updatePseudoData(
        indices[[m]], U[[inds[m, 1]]], U[[inds[m, 2]]],
        bias[[m]]$row$mu, bias[[m]]$col$mu
      )
      X[[m]][, 3] <- (
        xi - (1 / (1 + exp(-xi)) - origX[[m]]) / tau[m]
      )
    }

    for (m in which(likelihood == "poisson")) {
      xi <- p_updatePseudoData(
        indices[[m]], U[[inds[m, 1]]], U[[inds[m, 2]]],
        bias[[m]]$row$mu, bias[[m]]$col$mu
      )
      X[[m]][, 3] <- (
        xi - 1 / (1 + exp(-xi)) * (1 - origX[[m]] / log(1 + exp(xi))) / tau[m]
      )
    }

    par <- list(
      D = D, alpha = alpha, tau = tau,
      X = X, U = U, covU = covU, inds = inds, this = 1
    )

    if (!is.null(test)) {
      pred <- list()
      for (m in seq_len(M)) {
        indices_test <- matrix(as.integer(test[[m]][, 1:2]), ncol = 2)
        pred[[m]] <- p_updatePseudoData(
          indices_test, U[[inds[m, 1]]], U[[inds[m, 2]]],
          bias[[m]]$row$mu, bias[[m]]$col$mu
        )
      }
      for (m in which(likelihood == "bernoulli")) {
        pred[[m]] <- exp(pred[[m]])
        pred[[m]] <- pred[[m]] / (1 + pred[[m]])
      }
      for (m in which(likelihood == "poisson")) {
        pred[[m]] <- exp(pred[[m]])
        pred[[m]] <- log(1 + pred[[m]])
      }
      error <- rep(0, M)
      for (m in seq_len(M)) {
        error[m] <- sqrt(mean((test[[m]][, 3] - pred[[m]])^2))
      }
      errors[iter, ] <- error
      if ((opts$verbose > 0) && (iter %% 10 == 1)) {
        cat(paste0("Test error: ", paste(error, collapse = " "), "\n"))
      }
    }

    # Compute the cost function on training data
    if (opts$computeCost) {
      tcost <- 0
      for (m in seq_len(M)) {
        if (likelihood[m] == "gaussian") {
          # The likelihood term
          logtau <- digamma(a_tau[m]) - log(b_tau[m])
          tcost <- (
            tcost
            - dim(X[[m]])[1] * 0.5 * logtau
            + (b_tau[m] - a_tau[m]) * tau[m]
          )

          # KL divergence for tau
          temp <- (
            -lgamma(alpha_0t)
            + alpha_0t * log(beta_0t)
            + (alpha_0t - 1) * logtau
            - beta_0t * tau[m]
            + lgamma(a_tau[m])
            - a_tau[m] * log(b_tau[m])
            - (a_tau[m] - 1) * logtau
            + b_tau[m] * tau[m]
          )
          tcost <- tcost - temp
        } else {
          # Quadratic likelihood for non-Gaussian data
          temp <- sum(tau[m] * (X[[m]][, 3] - origX[[m]])^2 / 2)
          tcost <- tcost - temp
        }

        # Bias terms
        if (opts$useBias) {
          temp <- sum(
            -0.5 * log(bias[[m]]$row$nu)
            + 0.5 * log(bias[[m]]$row$scale)
            + (
              (bias[[m]]$row$nu + (bias[[m]]$row$m - bias[[m]]$row$mu)^2)
              / (2 * bias[[m]]$row$scale) - 0.5
            )
          )
          temp <- temp + sum(
            -0.5 * log(bias[[m]]$col$nu)
            + 0.5 * log(bias[[m]]$col$scale)
            + (
              (bias[[m]]$col$nu + (bias[[m]]$col$m - bias[[m]]$col$mu)^2)
              / (2 * bias[[m]]$col$scale) - 0.5
            )
          )
          temp <- (
            temp
            - 0.5 * log(bias[[m]]$row$scale)
            + 0.5 * log(1)
            + (
              (bias[[m]]$row$scale - (0 - bias[[m]]$row$m)^2)
              / (2 * 1^2) - 0.5
            )
          )
          temp <- (
            temp
            - 0.5 * log(bias[[m]]$col$scale)
            + 0.5 * log(1)
            + (
              (bias[[m]]$col$scale - (0 - bias[[m]]$col$m)^2)
              / (2 * 1^2) - 0.5
            )
          )
          tcost <- tcost + temp
        }
      }

      if (opts$method == "CMF") {
        for (i in seq_len(C)) {
          for (k in seq_len(K)) {
            logalpha <- digamma(a_ard[i]) - log(b_ard[i,k])

            # The U and covU terms
            temp <- (
              0.5 * dim(covU[[i]])[1] * logalpha
              - 0.5 * sum(covU[[i]][,k] + U[[i]][,k]^2) * alpha[i,k]
              + 0.5 * sum(log(covU[[i]][,k]))
            )
            tcost <- tcost - temp

            # The alpha part
            temp <- (
              -lgamma(alpha_0)
              + alpha_0 * log(beta_0)
              + (alpha_0 - 1) * logalpha
              - beta_0 * alpha[i, k]
              + lgamma(a_ard[i])
              - a_ard[i] * log(b_ard[i, k])
              - (a_ard[i] - 1) * logalpha + b_ard[i, k] * alpha[i, k]
            )
            tcost <- tcost - temp
          }
        }
      } else {
        for (k in seq_len(K)) {
          logalpha <- digamma(sum(a_ard)) - log(sum(b_ard))

          # The U and covU terms
          for (i in seq_len(C)) {
            temp <- (
              0.5 * dim(covU[[i]])[1] * logalpha
              - 0.5 * sum(covU[[i]][, k] + U[[i]][, k] ^ 2) * alpha[i, k]
              + 0.5 * sum(log(covU[[i]][, k]))
            )
            tcost <- tcost - temp
          }

          # The alpha part
          temp <- (
            -lgamma(alpha_0)
            + alpha_0 * log(beta_0)
            + (alpha_0 - 1) * logalpha
            - beta_0 * alpha[i, k]
            + lgamma(sum(a_ard))
            - sum(a_ard) * log(sum(b_ard[, k]))
            - (sum(a_ard) - 1) * logalpha
            + sum(b_ard[, k]) * alpha[i, k]
          )
          tcost <- tcost - temp
        }
      }
      cost <- c(cost, tcost)
    }

  } # the main loop of the algorithm ends

  Uall <- U[[1]]
  if (C > 1) {
    for (i in 2:C) {
      Uall <- rbind(Uall, U[[i]])
    }
  }

  # Reconstructed datasets
  if (!is.null(test)) {
    out <- list()
    for (m in seq_len(M)) {
      indices_test <- matrix(as.integer(test[[m]][, 1:2]), ncol = 2)
      xi <- p_updatePseudoData(
        indices_test, U[[inds[m, 1]]], U[[inds[m, 2]]],
        bias[[m]]$row$mu, bias[[m]]$col$mu
      )
      out[[m]] <- unname(cbind(indices_test, xi))
    }
    for (m in which(likelihood == "bernoulli")) {
      out[[m]][, 3] <- exp(out[[m]][, 3])
      out[[m]][, 3] <- out[[m]][, 3] / (1 + out[[m]][, 3])
    }
    for (m in which(likelihood == "poisson")) {
      out[[m]][, 3] <- exp(out[[m]][, 3])
      out[[m]][, 3] <- log(1 + out[[m]][, 3])
    }
  } else {
    out <- NULL
  }

  # return the output of the model as a list
  return(list(
    U = U,
    covU = covU,
    tau = tau,
    alpha = alpha,
    cost = cost,
    inds = inds,
    errors = errors,
    bias = bias,
    D = D,
    K = K,
    Uall = Uall,
    items = items,
    out = out,
    M = M,
    likelihood = likelihood,
    opts = opts
  ))
}

#' Internal function for computing the gradients
#'
#' @param r ?
#' @param par ?
#' @param stochastic Whether or not to perform updates on a subsample
#'
#' @return Gradient
p_gradUsparseWrapper <- function(r, par, stochastic = FALSE) {
  g <- matrix(r, nrow = par$D[par$this])
  cur <- g
  for (j in seq_len(par$D[par$this])) {
    g[j, ] <- g[j, ] * par$alpha[par$this, ]
  }

  for (m in which(par$inds[, 1] == par$this)) {
    v2 <- par$inds[m, 2]
    if (stochastic && m == 1) {
      part <- sample(nrow(par$X[[m]]), round(nrow(par$X[[m]]) / 10))
      # NOTE: p_gradUsparse directly modifies g
      p_gradUsparse(
        par$X[[m]][part, ], g, cur, par$U[[v2]], par$covU[[v2]], 1,
        par$tau[m], par$bias[[m]]$row$mu, par$bias[[m]]$col$mu
      )
    } else {
      p_gradUsparse(
        par$X[[m]], g, cur, par$U[[v2]], par$covU[[v2]], 1,
        par$tau[m], par$bias[[m]]$row$mu, par$bias[[m]]$col$mu
      )
    }
  }

  for (m in which(par$inds[, 2] == par$this)) {
    v1 <- par$inds[m, 1]
    if (stochastic && m == 1) {
      part <- sample(nrow(par$X[[m]]), round(nrow(par$X[[m]]) / 10))
      p_gradUsparse(
        par$X[[m]][part, ], g, cur, par$U[[v1]], par$covU[[v1]], 2,
        par$tau[m], par$bias[[m]]$row$mu, par$bias[[m]]$col$mu
      )
    } else {
      p_gradUsparse(
        par$X[[m]], g, cur, par$U[[v1]], par$covU[[v1]], 2,
        par$tau[m], par$bias[[m]]$row$mu, par$bias[[m]]$col$mu
      )
    }
  }

  return(g)
}

#' Internal function for checking whether the input is in the right format
#'
#' @param mat An input matrix of class `matrix`
#' @param max_row Maximum row index for `mat`
#' @param max_col Maximum column index for `mat`
#'
#' @return `TRUE` if the input is in coordinate/triplet format.
#'         `FALSE` otherwise.
p_check_sparsity <- function(mat, max_row, max_col) {
  if (is.matrix(mat)) {
    if (
      ncol(mat) != 3
      || min(mat[, 1:2]) < 1
      || max(mat[, 1]) > max_row
      || max(mat[, 2]) > max_col
    ) {
      cat("Matrix not in coordinate/triplet format")
      return(FALSE)
    } else {
      return(TRUE)
    }
  }

  cat("Input not of class `matrix`")
  return(FALSE)
}

#' Conversion from matrix to coordinate/triplet format
#'
#' The CMF code requires inputs to be speficied in a specific
#' sparse format. This function converts regular R matrices
#' into that format.
#'
#' The element `X[i, j]` on the `i`-th row and `j`-th column is represented
#' as a triple `(i, j, X[i,k])`. The input for CMF is then a matrix
#' where each row specifies one element, and hence the representation
#' is of size `N x 3`, where `N` is the total number of observed entries.
#'
#' In the original input matrix the missing entries should be marked
#' as `NA`. In the output they will be completely omitted.
#'
#' Even though this format reminds the representation often used
#' for representing sparse matrices, it is important to notice that
#' observed zeroes are retained in the representation. The
#' elements missing from this representation are considered unknown,
#' not zero.
#'
#' @param orig A matrix of class `matrix`
#' @return The input matrix in triplet/coordinate format.
#'
#' @author Arto Klami and Lauri Väre
#' @seealso [triplets_to_matrix()]
#' @examples
#'
#' x <- matrix(c(1, 2, NA, NA, 5, 6), nrow = 3)
#' triplet <- matrix_to_triplets(x)
#' print(triplet)
#'
#' @export
matrix_to_triplets <- function(orig) {
  triplets <- matrix(0, nrow = length(which(!is.na(orig))), ncol = 3)
  count <- 1

  for (y in seq_len(nrow(orig))) {
    for (x in seq_len(ncol(orig))) {
      if (!is.na(orig[y, x])) {
        triplets[count, ] <- c(y, x, orig[y, x])
        count <- count + 1
      }
    }
  }

  return(triplets)
}

#' Conversion from triplet/coordinate format to matrix
#'
#' This function is the inverse of [matrix_to_triplets()].
#' It converts a matrix represented as a set of triplets into
#' an object of the class `matrix`. The missing entries
#' (the ones not present in the triplet representation) are
#' filled in as `NA`.
#'
#' See [matrix_to_triplets()] for a description of the
#' representation.
#'
#' @param triplets A matrix in triplet/coordinate format
#' @return The input matrix as a normal matrix of class `matrix`
#' @author Arto Klami and Lauri Väre
#' @seealso [matrix_to_triplets()]
#' @examples
#'
#' x <- matrix(c(1, 2, NA, NA, 5, 6), nrow = 3)
#' triplet <- matrix_to_triplets(x)
#' print(triplet)
#' xnew <- triplets_to_matrix(triplet)
#' print(xnew)
#'
#' @export
triplets_to_matrix <- function(triplets) {
  mat <- matrix(NA, nrow = max(triplets[, 1]), ncol = max(triplets[, 2]))

  for (t in seq_len(nrow(triplets))) {
    mat[triplets[t, 1], triplets[t, 2]] <- triplets[t, 3]
  }

  return(mat)
}
