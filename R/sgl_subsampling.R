#
#     Description of this R script:
#     R interface to generic sparse group lasso subsampling procedure
#
#     Intended for use with R.
#     Copyright (C) 2014 Martin Vincent
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

# Tell R that 'task' (used in foreach) is global variable -- R check will complain if not
globalVariables('task')

#' Generic sparse group lasso subsampling procedure
#'
#' Support the use of multiple processors.
#'
#' @param module_name reference to objective specific C++ routines.
#' @param PACKAGE name of the calling package.
#' @param data a list of data objects -- will be parsed to the specified module.
#' @param parameterGrouping grouping of parameters, a vector of length \eqn{p}. Each element of the vector specifying the group of the parameters in the corresponding column of \eqn{\beta}.
#' @param groupWeights the group weights, a vector of length \code{length(unique(parameterGrouping))} (the number of groups).
#' @param parameterWeights a matrix of size \eqn{q \times p}.
#' @param alpha the \eqn{\alpha} value 0 for group lasso, 1 for lasso, between 0 and 1 gives a sparse group lasso penalty.
#' @param lambda the lambda sequence for the regularization path, a vector or a list of vectors (of the same length) with the lambda sequence for the subsamples.
#' @param training a list of training samples, each item of the list corresponding to a subsample.
#' Each item in the list must be a vector with the indices of the training samples for the corresponding subsample.
#' The length of the list must equal the length of the \code{test} list.
#' @param test a list of test samples, each item of the list corresponding to a subsample.
#' Each item in the list must be vector with the indices of the test samples for the corresponding subsample.
#' The length of the list must equal the length of the \code{training} list.
#' @param collapse if \code{TRUE} the results for each subsample will be collapse into one result (this is useful if the subsamples are not overlapping)
#' @param max.threads Deprecated (will be removed in 2018),
#' instead use \code{use_parallel = TRUE} and registre parallel backend (see package 'doParallel').
#' The maximal number of threads to be used.
#' @param use_parallel If \code{TRUE} the \code{foreach} loop will use \code{\%dopar\%}. The user must registre the parallel backend.
#' @param algorithm.config the algorithm configuration to be used.
#' @return
#' \item{responses}{content will depend on the C++ response class}
#' \item{features}{number of features used in the models}
#' \item{parameters}{number of parameters used in the models}
#' \item{lambda}{the lambda sequences used (a vector or list of length \code{length(training)}).}
#' @author Martin Vincent
#' @useDynLib sglOptim, .registration=TRUE
#' @importFrom utils packageVersion
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @export
sgl_subsampling <- function(module_name, PACKAGE,
	data,
	parameterGrouping,
	groupWeights,
	parameterWeights,
	alpha,
	lambda,
	training,
	test,
	collapse = FALSE,
	max.threads = NULL,
	use_parallel = FALSE,
	algorithm.config = sgl.standard.config) {

	# deprecated warnings
	if( ! is.null(max.threads) ) {
		warning("Argument 'max.threads' is deprecated (will be removed in 2018). \n    Alternative: Registre parallel backend (see package doParallel) and use 'use_parallel = TRUE' ")
	}

	# Check training and test consistency:
	if( ! is.list(training) | ! is.list(test)) {
		stop("The arguments traning and test should be lists")
	}

	if(length(training) != length(test)) {
		stop("The length of the lists traning and test should be equal")
	}

	# Check if collapse = TRUE is valid
	if(	collapse == TRUE
			&& (length(unique(unlist(test))) != data$n.samples || length(unlist(test)) != data$n.samples)) {
		stop("collapse = TRUE is not supported unless the test subsamples divide the samples into disjoint sets")
	}

	# Make lambda a list or validate
	if(is.list(lambda)) {
		if( ! all(sapply(lambda, function(x) length(x) == length(lambda[[1]])))) {
			stop("all lambda sequences must have same length")
		}

		if(length(lambda) != length(training)) {
			stop("The length of the lists lambda, traning and test should be equal")
		}

		lambda.list <- lambda

	} else {
		lambda.list <- replicate(length(training), lambda, simplify = FALSE)
	}

	training <- lapply(training, sort)
	test <- lapply(test, sort)

	call_sym <- paste(module_name, "sgl_subsampling", sep="_")

	# Registre parallel backend
	# This is only to make max.threads work -  remov in 2018
	if( ! is.null(max.threads) && max.threads > 1) {
		cl <- makeCluster(max.threads)
		registerDoParallel(cl)
	}

	if( use_parallel || ( ! is.null(max.threads) && max.threads != 1 ) ) {

		rawres <- foreach(task=1:length(training),
			.packages = PACKAGE) %dopar% {

           test_data <- subsample(data, test[[task]])
			train_data <- subsample(data, training[[task]])

			# Prapare arguments
			args <- prepare.args(
				data = train_data,
				parameterGrouping = parameterGrouping,
				groupWeights = groupWeights,
				parameterWeights = parameterWeights,
				alpha = alpha,
				test_data = test_data)

			.Call(call_sym, PACKAGE = PACKAGE,
				args$data,
				args$test_data,
				args$block.dim,
				args$groupWeights,
				args$parameterWeights,
				args$alpha,
				lambda.list[[task]],
				algorithm.config)
			}

	} else {

		rawres <- foreach(task=1:length(training)) %do% {

       test_data <- subsample(data, test[[task]])
	    train_data <- subsample(data, training[[task]])

		# Prapare arguments
		args <- prepare.args(
			data = train_data,
			parameterGrouping = parameterGrouping,
			groupWeights = groupWeights,
			parameterWeights = parameterWeights,
			alpha = alpha,
			test_data = test_data)

		.Call(call_sym, PACKAGE = PACKAGE,
				args$data,
				args$test_data,
				args$block.dim,
				args$groupWeights,
				args$parameterWeights,
				args$alpha,
				lambda.list[[task]],
				algorithm.config)
			}
	}

	if( ! is.null(max.threads) && max.threads > 1) {
 		stopCluster(cl)
	}

	# formating responses
	res <- list()

	res$responses <- .format_responses(
		rawres,
		collapse,
		length(training),
		length(lambda))

	res$features <- t(sapply(rawres, function(x) x$features))
	res$parameters <- t(sapply(rawres, function(x) x$parameters))

	# Sample names
	sample.names <- data$sample.names

	if(collapse == TRUE) {

		# Reorder responses and set sample names
		sample.order <- order(unlist(test))
		res$responses <- lapply(res$responses, function(x) .order_response(x, sample.order))
		res$responses <- lapply(res$responses, function(x) .set_sample_names(x, sample.names))

	} else {
		#Set sample names
		res$responses <- lapply(res$responses, function(x) lapply(1:length(x), function(i) .set_sample_names(x[[i]], sample.names[test[[i]]])))
	}

	# Names
	rownames(res$features) <- paste("subsample", 1:length(training))
	rownames(res$parameters) <- paste("subsample", 1:length(training))

	res$lambda <- lambda.list

	# Set version, type and class and return
	res$sglOptim_version <- packageVersion("sglOptim")
	res$type <- "subsampling"
	class(res) <- "sgl"

	return(res)
}

.format_responses <- function(x, collapse, n_subsamples, n_lambda) {
	if( ! collapse) {

		r <- lapply(1:length(x[[1]]$responses), function(k)

			if(is.list(x[[1]]$responses[[k]])) {

			 return( lapply(1:n_subsamples,
					function(i) lapply(1:n_lambda,
						function(j) x[[i]]$responses[[k]][[j]])))

			} else {

				return( lapply(1:n_subsamples,
							function(j) x[[j]]$responses[[k]]) )
			})

	} else {

		r <- lapply(1:length(x[[1]]$responses), function(k)

			if (is.list(x[[1]]$responses[[k]])) {

			 return( lapply(1:n_lambda,
					function(i) .collapse_responses(lapply(1:n_subsamples,
						function(j) x[[j]]$responses[[k]][[i]]))) )

			} else {

				return( .collapse_responses(lapply(1:n_subsamples,
							function(j) x[[j]]$responses[[k]])) )

			})
	}

	names(r) <- names(x[[1]]$responses)

	return(r)
}

.collapse_responses <- function(responses) {
	if(is.matrix(responses[[1]])) {
		return( do.call(cbind, responses) )
	} else {
		stop("unknown response class")
	}
}

.set_sample_names <- function(response, sample.names) {

	if(is.list(response) && is.list(response[[1]]) && ! is.null(response[[1]]$parameters) && response[[1]]$parameters == TRUE) {
		names(response) <- sample.names
		return(response)
	}

	if(is.list(response)) {
		return(lapply(response, function(x) .set_sample_names(x, sample.names)))
	}

	if(is.matrix(response)) {
		colnames(response) <- sample.names
		return(response)
	}

	if(is.vector(response)) {
		names(response) <- sample.names
		return(response)
	}

	stop("Unknown response class")

}

.order_response <- function(response, sample.order) {

	if(is.list(response) && is.list(response[[1]]) && ! is.null(response[[1]]$parameters) && response[[1]]$parameters == TRUE) {
		return(response[sample.order])
	}

	if(is.list(response)) {
		return(lapply(response, function(x) .order_response(x, sample.order)))
	}

	if(is.matrix(response)) {
		return(response[,sample.order, drop = FALSE])
	}

	if(is.vector(response)) {
		return(response[sample.order])
	}

	stop("Unknown response class")
}
