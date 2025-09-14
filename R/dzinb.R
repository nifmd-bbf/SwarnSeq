#' @title Probability Density for the Zero-Inflated Negative Binomial Distribution
#' @description
#' This function calculates the probability density for a Zero-Inflated Negative Binomial (ZINB) distribution.
#' @param x This is a vector of non-negative counts for a gene across all the cells.
#' @param size This is the over-dispersion parameter of the negative binomial distribution.
#' @param mu The mean parameter of the negative binomial distribution.
#' @param rho The zero-inflation parameter, a probability value between 0 and 1.
#' @param log A logical value. If TRUE, this function returns the logarithm of the density. Otherwise, it returns the density itself.
#' @returns
#' The probability density.
#' @export
#' @examples
#' dzinb_results <- SwarnSeq::dzinb(x = c(1, 2, 5, 4, 7, 0), size = 1, mu = 3, rho = .4, log = TRUE)
dzinb <- function(x, size, mu, rho, log) {
    show.custom.warning <- function(x) warning(x)
    if (sum(is.na(x)) > 0) {
        show.custom.warning("NA detected in x..")
        return(invisible(NULL))
    }
    if (sum(x < 0)) {
        show.custom.warning("Negative values detected in x..")
        return(invisible(NULL))
    }
    if (size < 0 || sum(mu < 0) > 0 || sum(rho < 0) > 0) {
        show.custom.warning("All parameters should be positive..")
        return(invisible(NULL))
    }
    if (sum(rho < 0) > 0 || sum(rho > 1) > 0) {
        show.custom.warning("The probability of structural zeros must be [0, 1].")
        return(invisible(NULL))
    }
    if (!is.logical(log) || length(log) != 1) {
        show.custom.warning("Bad input for argument 'log'..")
        return(invisible(NULL))
    }
    out <- rho * (x == 0) + (1 - rho) * dnbinom(x, size = size, mu = mu, log = FALSE)
    if (log == TRUE) {
        return(log(out))
    }else {
        return(out)
    }
}
