#' Simulated TV-VAR dataset
#'
#' Artificial dataset generated from a TV-VAR(1) model with 2 variables
#' and 1 factor. The parameters used in the simulation are:
#'
#' \describe{
#'   \item{dim.VAR}{2}
#'   \item{lag.order}{1}
#'   \item{number.factors}{1}
#'   \item{T}{1000 observations}
#'   \item{A2}{0.1}
#'   \item{B2}{0.75}
#'   \item{pphi}{0.95}
#'   \item{omega}{\eqn{\begin{bmatrix}0.3 & 0.2 \\ 0.2 & 0.3\end{bmatrix}}}
#'   \item{Phi.c}{\eqn{\begin{bmatrix}0.3 & 0.15 \\ 0.0 & 0.3\end{bmatrix}}}
#'   \item{Phi.f}{\eqn{\begin{bmatrix}0.2 & 0.0 \\ 0.0 & 0.2\end{bmatrix}}}
#' }
#'
#' @format A list with:
#' \describe{
#'   \item{data}{T x N matrix of simulated observations.}
#'   \item{Phi.f.array}{Array of factor loadings used.}
#'   \item{u}{Innovations.}
#'   \item{factors}{Latent factor series.}
#' }
#'
"simdata"

#' Simulated Sparse tvVAR(2) Dataset
#'
#' Artificial dataset generated from a sparse time-varying VAR(2) model
#' with 2 observed variables and 2 latent factors.  
#' This dataset illustrates how sparsity in the factor loadings \eqn{\Phi_f}
#' affects estimation and model selection with penalization.
#'
#' @details
#' Simulation setup:
#' \describe{
#'   \item{N}{2 observed variables.}
#'   \item{p}{2 lags.}
#'   \item{r}{2 latent factors.}
#'   \item{A2, B2}{Scalar BEKK(1,1) parameters controlling volatility dynamics
#'     (\eqn{A^2 = 0.1}, \eqn{B^2 = 0.75}).}
#'   \item{pphi}{Diagonal factor AR(1) coefficients (\eqn{0.95}).}
#'   \item{omega}{Innovation covariance matrix
#'     \eqn{\begin{bmatrix}0.3 & 0.2 \\ 0.2 & 0.3\end{bmatrix}}.}
#'   \item{Phi.c}{Static VAR coefficients
#'     \eqn{\begin{bmatrix}0.30 & 0.10 & 0 & 0.10 \\
#'                          0 & 0.25 & 0.10 & 0\end{bmatrix}}.}
#'   \item{Phi.f}{Time-varying (factor-driven) loadings across lags and factors:
#'     a 3D array \eqn{[N \times (N p) \times r]} with sparsity pattern
#'      \eqn{: \Phi_f^{(1)} = \begin{bmatrix}
#'         0 & 0.25 & 0 & 0 \\
#'         0 & 0 & 0 & -0.20
#'       \end{bmatrix},
#'       \Phi_f^{(2)} =
#'       \begin{bmatrix}
#'         0 & 0 & 0 & 0.15 \\
#'         0 & 0 & 0 & 0
#'       \end{bmatrix}.
#'     }
#'   }
#'  \item{burn_in}{500 initial observations discarded.}
#' }
#' @format
#' A list with the following elements:
#' \describe{
#'   \item{Y}{Numeric matrix (\eqn{T \times N}) of simulated observations.}
#'   \item{factors}{Numeric matrix (\eqn{T \times r}) of latent factors.}
#'   \item{shocks}{Numeric matrix (\eqn{T \times N}) of structural shocks.}
#' }
#'
#' @source Simulated via \code{\link{tvsim}}.
#'
#' @examples
#' data(simsparse)
#' str(simsparse)
#' matplot(simsparse$Y, type = "l", col = 1:2, lty = 1,
#'         main = "Simulated Sparse tvVAR(2) Data",
#'         ylab = "Series values", xlab = "Time")
"simsparse"