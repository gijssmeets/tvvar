#' Penalized estimation (ECM) with optional adaptive weights
#'
#' @param data            numeric matrix (T x N)
#' @param p               integer VAR lag order
#' @param r               integer number of factors
#' @param zero.mean       logical; if TRUE, intercepts fixed at 0
#' @param lambda_penalty  numeric scalar; L1 penalty level for Phi.f
#' @param penalty_type    "adaptive" or "regular"
#' @param Phi.f.structure optional free-pattern for Phi^f (3D [N,N,r], list of r N x N, or N x (N*r) matrix)
#' @param init            "default", "random", or "custom"
#' @param init_list       named list for custom init (A, B, phi_r, Omega, Phi_c, Phi_f)
#'
#' @return tvfit list (same fields as unpenalized_estimate output)
#' @export
penalized_estimate <- function(data,
                               p = 1,
                               r = 1,
                               zero.mean = TRUE,
                               lambda_penalty = 0.01,
                               penalty_type = c("adaptive", "regular"),
                               Phi.f.structure = NULL,
                               init = c("default", "random", "custom"),
                               init_list = NULL) {
  stop("penalized_estimate() is disabled in this lite build (no 'apg'). ",
       "Install the full tvvar package to enable penalized estimation.")
}