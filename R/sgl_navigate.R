#
#     Description of this R script:
#     Routines for navigating sgl objects
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

# S3 functions:

#' @title Generic function for computing error rates
#'
#' @description
#' Compute and returns an error rate for each model contained in \code{x}.
#' See details for generic use cases.
#'
#' @details
#' The following generic use case should be supported (see for example \pkg{msgl} package for an implementation):
#'
#' \enumerate{
#' \item With \code{fit} a sgl fit object with models estimated using \code{x} data, the code
#'
#' \code{Err(fit, x)}
#'
#' should return a vector with the \emph{training errors} of the models.
#'
#' \item With \code{x.new} a new data set with known responses \code{response.new}, the code
#'
#' \code{Err(fit, x.new, response.new)}
#'
#' should return a vector with the errors of the models when applied to the new data set.
#
#' \item With \code{fit.cv} a sgl cross validation object, the code
#'
#' \code{Err(fit.cv)}
#'
#' should return a vector with estimates of the \emph{expected generalization errors} of the models (i.e. the cross validation errors).
#'
#' \item If subsampling is supported then, with \code{fit.sub} a sgl subsampling object, the code
#'
#' \code{Err(fit.sub)}
#'
#' should return a matrix with the test errors (each column corresponding to a model, i.e. rows corresponds to tests).
#' }
#'
#' @seealso compute_error
#' @param object an object
#' @param data a data object
#' @param response a response object
#' @param ... additional parameters (optional)
#' @return
#' a vector of length \code{nmod(object)} or a matrix with \code{nmod(object)} columns containing error rates for the models
#' @author Martin Vincent
#' @export
Err <- function(object, data, response, ... ) UseMethod("Err")


#' @title Helper function for computing error rates
#' @description
#' This function can be used to compute error rates.
#' It is consist with the use cases of the \code{Err} genetic function.
#' (see \pkg{msgl} package for an example of how to use this function)
#' @param object a object containing responses
#' @param data a data object
#' @param response.name the name of the response, if \code{response.name = NULL} then \code{object} will be treated as the response.
#' @param response the response
#' @param loss the loss function
#' @return a vector with the computed error rates
#'
#' @author Martin Vincent
#' @importFrom stats predict
#' @export
compute_error <- function(object, data = NULL, response.name, response, loss) {

	if(!is.null(data)) {
		return(compute_error(object = predict(object, data), data = NULL, response.name = response.name, response = response, loss = loss))
	}

	if(is.null(response.name) || any(names(object) == response.name)) {

		if(is.null(response.name)) {
			r <- object
		} else {
			r <- object[[response.name]]
		}

		if(is.list(r) && is.list(response)) {
			return(t(sapply(1:length(r), function(i) compute_error(object = r[[i]], data = NULL, response.name = NULL, response = response[[i]], loss = loss))))

		} else if(is.list(r)) {
			return(sapply(r, function(m) loss(m, response)))

		} else if(is.matrix(r)) {
			return(apply(r, 2, FUN = function(v) loss(v, response)))

		} else if(is.vector(r))  {
			return(loss(r, response))

		} else {
			stop("Unknown response type")
		}
	}

	stop(paste("response '", response.name, "' not found", sep =""))

}

#' @title Extracts nonzero features
#'
#' @description
#' Generic function for extracting nonzero features.
#'
#' @param object an object
#' @param ... additional parameters (optional)
#' @return a list of length \code{nmod(x)} containing the nonzero features of the models.
#'
#' @author Martin Vincent
#' @export
features <- function(object, ...) UseMethod("features")

#' @title Extract feature statistics
#' @description
#' Generic function for extracting feature statistics.
#'
#' @param object an object
#' @param ... additional parameters (optional)
#' @return an object conting the computed statistics.
#'
#' @author Martin Vincent
#' @export
features_stat <- function(object, ...) UseMethod("features_stat")

#' @title Extracts nonzero parameters
#'
#' @description
#' Generic function for extracting nonzero parameters for each model.
#'
#' @param object an object
#' @param ... additional parameters (optional)
#' @return a list of length \code{nmod(x)} containing the nonzero parameters of the models.
#'
#' @author Martin Vincent
#' @export
parameters <- function(object, ...) UseMethod("parameters")

#' @title Extract parameter statistics
#' @description
#' Generic function for extracting parameter statistics.
#'
#' @param object an object
#' @param ... additional parameters (optional)
#' @return an object conting the computed statistics.
#'
#' @author Martin Vincent
#' @export
parameters_stat <- function(object, ...) UseMethod("parameters_stat")

#' @title Number of models used for fitting
#'
#' @description
#' Generic function for counting the number of models used for fitting the object.
#' Returns the number of models used for fitting.
#' However, note that the objects returned by \code{msgl.cv} and {msgl.subsampling} does not contain any models even though \code{nmod} returns a nonzero number.
#'
#' @param object an object
#' @param ... additional parameters (optional)
#' @return the number of models used when fitting the object \code{x}.
#'
#' @author Martin Vincent
#' @export
nmod <- function(object, ...) UseMethod("nmod")

#' @title Index of best model
#'
#' @description
#' Returns the index of the best model
#' @param object an sgl object
#' @param ... additional parameters (optional)
#' @return index of the best model.
#'
#' @author Martin Vincent
#' @export
best_model <- function(object, ...) UseMethod("best_model")

#' @title Extract fitted models
#'
#' @description
#' Generic function for extracting the fitted models.
#' Returns the fitted models.
#'
#' @param object an object
#' @param index a vector of indices of the models to be returned
#' @param ... additional parameters (optional)
#' @return a list of length \code{length(index)} containing the models
#'
#' @author Martin Vincent
#' @export
models <- function(object, index, ...) UseMethod("models")


#' @title Extracting nonzero features
#' @description
#' Extract the nonzero features of the fitted models
#'
#' @param object a sgl object
#' @param ... ignored
#' @return a list of vectors containing the nonzero features (that is nonzero columns of the \eqn{beta} matrices)
#'
#' @author Martin Vincent
#' @export
features.sgl <- function(object, ...) {

	if(is.null(object$beta)) {
		stop("object contains no models")
	}

	if(is.null(colnames(object$beta[[1]])) || any(duplicated(colnames(object$beta[[1]])))) {
		res <- lapply(object$beta, function(beta) which(colSums(beta != 0) != 0))
	} else {
		res <- lapply(object$beta, function(beta) colnames(beta)[colSums(beta != 0) != 0])
	}

	return(res)
}

#' @title Extract feature statistics
#'
#' @description
#' Extracts the number of nonzero features (or group) in each model.
#'
#'
#' @param object an object
#' @param ... ignored
#' @return a vector of length \code{nmod(x)} or a matrix containing the number of nonzero features (or group) of the models.
#'
#' @author Martin Vincent
#' @export
features_stat.sgl <- function(object, ...) {

	if(object$type == "fit") {
		return(sapply(features(object), length))
	} else if (object$type == "cv") {
		return(object$features)
	} else if (object$type == "subsampling") {
		return(object$features)
	} else {
		stop("Unknown sgl object type")
	}

}

#' @title Extracting nonzero parameters
#' @description
#' Extract the nonzero parameters in each model.
#' Only the parameters of nonzero features (columns of the \eqn{beta} matrices) are returned.
#'
#' @param object a sgl object
#' @param ... ignored
#' @return a list of vectors containing the nonzero parameters (that is nonzero entries of the \eqn{beta} matrices)
#'
#' @author Martin Vincent
#' @export
parameters.sgl <- function(object, ...) {

	if(is.null(object$beta)) {
		stop("object contains no models")
	}

	tmp <- features(object)
	res <- sapply(1:length(object$beta), function(i) object$beta[[i]][,tmp[[i]], drop = FALSE] != 0)

	return(res)
}

#' @title Extracting parameter statistics
#'
#' @description
#' Extracts the number of nonzero parameters in each model.
#'
#' @param object an object
#' @param ... ignored
#' @return a vector of length \code{nmod(x)} or a matrix containing the number of nonzero parameters of the models.
#'
#' @author Martin Vincent
#' @export
parameters_stat.sgl <- function(object, ...) {

	if(object$type == "fit") {
		return(sapply(parameters(object), length))
	} else if (object$type == "cv") {
		return(object$parameters)
	} else if (object$type == "subsampling") {
		return(object$parameters)
	} else {
		stop("Unknown sgl object type")
	}
}

#' @title Returns the number of models in a sgl object
#' @description
#' Returns the number of models used for fitting.
#' @details
#' Note that cv and subsampling objects does not conating any models even though nmod returns a positiv number.
#'
#' @param object a sgl object
#' @param ... ignored
#' @return the number of models in \code{object}
#'
#' @author Martin Vincent
#' @export
nmod.sgl <- function(object, ...) {
	if(is.list(object$lambda)) {
			return(length(object$lambda[[1]]))
	}

	return(length(object$lambda))
}

#' @title Index of best model
#'
#' @description
#' Returns the index of the best model, in terms of lowest error rate
#' @param object a sgl object
#' @param pkg name of calling package
#' @param ... additional parameters (optional)
#' @return index of the best model.
#'
#' @author Martin Vincent
#' @export
best_model.sgl <- function(object, pkg, ...) {

#NOTE should be consistent with sgl_print

	class(object) <- pkg

	if(object$type == "fit") {
    stop("not applicable to fit objects")

  } else if(object$type == "cv") {

    err <- Err(object)
    return(which(err == min(err))[1])

  } else if(object$type == "subsampling") {

    err <- Err(object)

    model.sel <- apply(err, 1, function(y) which(min(y) == y)[1])
    return(model.sel)
  }

  stop("Unknown type of sgl object")
}


#' @title Extract the estimated models
#'
#' @description
#' This function returns the estimated models (that is the \eqn{beta} matrices)
#'
#' @param object a sgl object
#' @param index indices of the models to be returned
#' @param ... ignored
#' @return a list of sparse matrices
#'
#' @author Martin Vincent
#' @export
models.sgl <- function(object, index = 1:nmod(object), ...) {

	if(is.null(object$beta)) {
		stop("object contains no models")
	}

	return(object$beta[index])
}

#' @title Extracting the nonzero coefficients
#' @description
#' This function returns the nonzero coefficients (that is the nonzero entries of the \eqn{beta} matrices)
#'
#' @param object a sgl object
#' @param index indices of the models
#' @param parameter name of the parameter (deafult is 'beta')
#' @param ... ignored
#' @return a list of with nonzero coefficients of the models
#'
#' @author Martin Vincent
#' @export
coef.sgl <- function(object, index = 1:nmod(object), parameter = "beta", ...) {

	if(is.null(object[[parameter]])) {
		stop("object contains no models")
	}

	return(get_coef(object[[parameter]], index))
}

#' @title Get the nonzero coefficients
#' @description
#' Extracting nonzero coefficients from list (of lists) of matrices
#'
#' @param object a list of lists of matrices or a list of matrices
#' @param index indices to be extracted from
#' @return a list (of lists) with the nonzero coefficients
#'
#' @author Martin Vincent
#' @export
get_coef <- function(object, index = 1:length(object)) {
    lapply(object[index], function(x) {
      if(is.list(x)) {
        return(get_coef(x))
      } else {
        return(x[,colSums(x != 0) != 0, drop = FALSE])
      }})
  }

#' @title Print information about sgl object
#' @description
#' Prints information about sgl object
#'
#' @param x a object of sgl family class
#'
#' @author Martin Vincent
#' @export
sgl_print <- function(x) {

	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
			"\n", sep = "")

	if(x$type == "fit") {

		cat("\nModels:\n\n")
		sel <- 1:5*floor(nmod(x)/5)

		feat <- sapply(features(x), length)
		para <- sapply(parameters(x), sum)

		print(data.frame(
						'Index: ' = sel,
						'Lambda: ' = x$lambda[sel],
						'Features: ' = print_with_metric_prefix(feat[sel]),
						'Parameters: ' = print_with_metric_prefix(para[sel]), check.names = FALSE),
				row.names = FALSE, digits = 2, right = TRUE)

		cat("\n")

	} else if(x$type == "cv") {

		cat("\nModels:\n\n")
		sel <- 1:5*floor(nmod(x)/5)

		err <- Err(x)
		feat <- colMeans(x$features)
		para <- colMeans(x$parameters)

		print(data.frame('Index: ' = sel,
						'Lambda: ' = x$lambda[sel],
						'Features: ' = print_with_metric_prefix(feat[sel]),
						'Parameters: ' = print_with_metric_prefix(para[sel]),
						'Error: ' = err[sel], check.names = FALSE),
				row.names = FALSE, digits = 2, right = TRUE)

		cat("\nBest model:\n\n")

		sel <- which(err == min(err))[1]

		print(data.frame('Index: ' = sel,
						'Lambda: ' = x$lambda[sel],
						'Features: ' = print_with_metric_prefix(feat[sel]),
						'Parameters: ' = print_with_metric_prefix(para[sel]),
						'Error: ' = err[sel], check.names = FALSE),
				row.names = FALSE, digits = 2, right = TRUE)

		cat("\n")

	} else if(x$type == "subsampling") {

		cat("\nBest models:\n\n")

		err <- Err(x)

		if(nrow(err) <= 5) {
			sel <- 1:nrow(err)
		} else {
			sel <- 1:5*floor(nrow(err)/5)
		}

		model.sel <- apply(err[sel,], 1, function(y) which(min(y) == y)[1])
		feat <- sapply(1:length(sel), function(i) x$features[sel[i], model.sel[i]])
		para <- sapply(1:length(sel), function(i) x$parameters[sel[i], model.sel[i]])

		print(data.frame('Subsample: ' = sel,
						'Model index: ' = model.sel,
						'Lambda: ' = sapply(1:length(sel), function(i) x$lambda[[ sel[i] ]][ model.sel[i] ]),
						'Features: ' = print_with_metric_prefix(feat),
						'Parameters: ' = print_with_metric_prefix(para),
						'Error: ' = sapply(1:length(sel), function(i) err[sel[i], model.sel[i]]), check.names = FALSE),
				row.names = FALSE, digits = 2, right = TRUE)

		cat("\n")

	} else if(x$type == "predict") {

		cat("\nPredictions for ")
		cat(x$n.samples)
		cat(" samples at ")
		cat(length(x$lambda))
		cat(" lambda values.\n\n")

		cat("Lambda values: ")
		if(length(x$lambda) >= 4) {
			cat(paste(c(round(x$lambda[1:3],3), "...", x$lambda[length(x$lambda)]), collapse = " "))
		} else {
			cat(paste(x$lambda, collapse = " "))
		}
		cat("\n")

		if( ! is.null(x$sample.names)) {
			cat("\nSamples: ")
			if(length(x$sample.names) >= 4) {
				cat(paste(c(x$sample.names[1:3], "...", x$sample.names[length(x$sample.names)]), collapse = " "))
			} else {
				cat(paste(x$sample.names, collapse = " "))
			}
			cat("\n\n")
		}

	} else {
		stop("Unknown type of sgl object")
	}
}


#' Print a numeric with metric prefix
#'
#' @param x numeric to be printed
#' @param digits number of significant digits
#'
#' @return a string
#'
#' @author Martin Vincent
#' @export
print_with_metric_prefix <- function(x, digits = 3) {

	if(length(x) > 1) {
		return(sapply(x, function(y) print_with_metric_prefix(y, digits = digits)))
	}

	metric_factor <-  c(1e+00, 1e+03, 1e+06, 1e+09)

	prefix <- c("", "k", "M", "G")

	sel <- max(which(x >= metric_factor))

	txt <- paste(round(x/metric_factor[sel], digits = digits), prefix[sel], sep="")

	return(txt)
}
