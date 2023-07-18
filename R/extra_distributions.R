## Extra distributions for extracting estimates of the EVZINB and EVINB
nbinomdist2 <- function (size = 10, prob, mu){
  if (missing(mu) & !missing(prob)) {
    if (!is.numeric(prob) || !is.numeric(size))
      stop("Parameters must be a numeric")
    if (prob < 0 || prob > 1)
      stop("prob must be in [0,1].")
    if (size != floor(size))
      stop("size must be an integer.")
    # x <- list(parameters = list(size = size, prob = prob),
    #           type = "Negative Binomial", support = list(from = 0,
    #                                                      to = Inf, by = 1))
  }
  else {
    if (!missing(mu) & missing(prob)) {
      if (!is.numeric(mu) || !is.numeric(size))
        stop("Parameters must be a numeric")
      # if (size != floor(size))
      #   stop("size must be an integer.")
      # x <- list(parameters = list(size = size, mu = mu),
      #           type = "Negative Binomial", support = list(from = 0,
      #                                                      to = Inf, by = 1))
    }
    else {
      stop("Either prob or mu has to be set.")
    }
  }
  mistr::new_dist(name = 'Negbin2',
           from = 0, to = Inf, by = 1,
           class = c("nbinomdist", "discrdist", "standist", "univdist",
                     "dist"))

}
