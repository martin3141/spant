#' Generate a gaussian pulse shape
#' @param angle pulse angle in degrees
#' @param n number of points to generate
#' @param trunc percentage truncation factor
#' @export
get_guassian_pulse <- function(angle, n, trunc = 1) {
  T = 1
  t <- seq(from = -T / 2, to = T / 2, length.out = n)
  beta <- -log(trunc / 100)
  B <- exp(-beta * (2 * t / T) ^ 2)
  B <- B / sum(B) * pi * angle / 180
  B
}