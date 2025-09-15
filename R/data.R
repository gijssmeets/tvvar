#' Macro dataset: IP, CPI, bond spread
#'
#' Monthly US data used in tvvar examples.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{date}{Date (character or Date)}
#'   \item{ip}{Industrial production index}
#'   \item{cpi}{Consumer Price Index}
#'   \item{bond_spread}{Bond yield spread}
#'   \item{recession}{Recession dummy}
#' }
"macro_data"

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