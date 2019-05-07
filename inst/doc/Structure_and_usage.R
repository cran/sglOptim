## ---- eval=FALSE---------------------------------------------------------
#  #' C interface
#  #'
#  #' @keywords internal
#  #' @export
#  msgl_dense_sgl_fit_R <- function(
#    data,
#    block_dim,
#    groupWeights,
#    parameterWeights,
#    alpha,
#    lambda,
#    idx,
#    algorithm.config) {
#  
#    .Call(msgl_dense_sgl_fit, PACKAGE = "msgl",
#          data,
#          block_dim,
#          groupWeights,
#          parameterWeights,
#          alpha,
#          lambda,
#          idx,
#          algorithm.config
#    )
#  }

