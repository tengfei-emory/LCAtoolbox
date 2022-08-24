#' Simulate a dataset which can be analyzed by SLTCA
#'
#' @description Simulate a dataset with longitudinal observations.
#' @param n Sample size.
#' @author Teng Fei. Email: tfei@emory.edu
#' @return Returns a data frame with 6 longitudinal features y.1 - y.6, including count (y.1 and y.2), binary (y.3 and y.4) and continuous (y.5 and y.6) type. Variable baselinecov is the baseline risk factor of latent classes. Variable latent is the true latent class labels.
#' @examples
#'
#' dat <- simulation.SLTCA(500)
#' @export


simulation.SLTCA <- function(n){
  output <- simulate(n)
  return(output)
}