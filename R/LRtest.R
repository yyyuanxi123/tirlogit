#' Likelihood ratio test for `tirlogit` models
#'
#' The `LRtest` function performs a likelihood ratio tests
#' for two nested `tirlogit` models, to check whether
#' one model provides a significantly
#' better fit to the data than the other model.
#'
#' @name LRtest
#' @param model1 a fitted model of class `tirlogit`.
#' @param model2 a fitted model of class `tirlogit`.
#' @return Summary of the likelihood ratio test, showing degrees of freedom,
#' log-likelihood values, chi-square statistic, and p-value.
#' @examples
#' ### Conduct a likelihood ratio test for nested models ###
#' data(mode_long, package = "tirlogit")
#' model1 <- tirlogit(rank ~ cost, data = mode_long, intercept = FALSE, reflevel = "car",
#'                    type = "rank", ties = TRUE, shape = "long", alts.names = "alt",
#'                    indexes = c("case", "taskid"))
#' model2 <- tirlogit(rank ~ cost + freq, data = mode_long, intercept = FALSE, reflevel = "car",
#'                    type = "rank", ties = TRUE, shape = "long", alts.names = "alt",
#'                    indexes = c("case", "taskid"))
#' LRtest(model1, model2)
#' @export
LRtest <- function(model1, model2){
  df_model1 <- length(model1$coefficients)
  df_model2 <- length(model2$coefficients)
  loglike_model1 <- model1$LogLikelihood
  loglike_model2 <- model2$LogLikelihood
  df_21 <- df_model2 - df_model1
  if (df_21 > 0){
    df <- df_21
    chisq <- -2*(loglike_model1 - loglike_model2)
  }else{
    df <- -df_21
    chisq <- -2*(loglike_model2 - loglike_model1)
  }

  p_value <- 1 - pchisq(chisq, df)
  chisq <- round(chisq, 6)
  p_value <- round(p_value, 6)
  LR_results <- data.frame(
    Df = c(df_model1, df_model2),
    LogLik = c(loglike_model1, loglike_model2),
    Df2 = c("", abs(df_21)),
    Chisq = c("", chisq),
    PrChisq = c("", p_value)
  )

  colnames(LR_results) <- c("#Df", "LogLik", "Df", "Chisq", "Pr(>Chisq)")

  cat("Likelihood ratio test\n\n")
  cat("Model 1:")
  print(model1$call$formula)
  cat("Model 2:")
  print(model2$call$formula)
  print(LR_results)
  if(p_value <=0.1){
    cat("---\n")
    cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")}
}

