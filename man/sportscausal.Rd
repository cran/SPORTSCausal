\name{sportscausal}
\alias{sportscausal}
\title{Time series causal inference of Randomized Controlled Trial (RCT) under spillover effect}
\description{
'SPORTSCausal' produces treatment effect and spillover effect estimation from responses of experiment group and control group.
}
\usage{
sportscausal(y.exp, y.con, pre.period, post.period, is.plot = TRUE, 
  model.select = "AIC", max.p = 3, max.d = 3, max.q = 3, feature = NULL)
}
\arguments{
  \item{y.exp}{Response of experiment group, from pre-treatment to post-treatment}
  \item{y.con}{Response of control group, from pre-treatment to post-treatment}
  \item{pre.period}{Time period before the treatment}
  \item{post.period}{Time period during the treatment}
  \item{is.plot}{If \code{is.plot = TRUE}, by default, a pdf containing summary figures will be returned to the current working directory as \code{getwd()}}
  \item{model.select}{Model used to predict the time series without treatment. If \code{model.select = "AIC"}, by default, the ARIMA model using AIC selection would be applied. If \code{model.select = "CV"}, the ARIMA model using cross validation would be applied. If \code{model.select = "lstm"}, the Long Short-Term Memory model would be applied}
  \item{max.p}{The max number of autoregressive terms in ARIMA model, by default \code{max.p = 3}}
  \item{max.d}{The max number of nonseasonal differences needed for stationarity in ARIMA model, by default \code{max.d = 3}}
  \item{max.q}{The max number of lagged forecast errors in the prediction equation in ARIMA model, by default \code{max.p = 3}}
  \item{feature}{The covariate matrix associated with the response. By default, \code{feature = NULL} but can be non-null when \code{model.select = "lstm"}}
}
\details{In the presense of spillover effect, the response of control group could be interferenced by the treatment. In order to seprate the treatment effect and spillover effect, \code{sportscausal} uses ARIMA model or LSTM model to predict the response behavior without treatment. The point estimator and significance of both effect follow using Bayesian Structrual Time Series (BSTS) model.}
\value{
  \item{est.treatment}{Information of treatment effect estimation, containing point estimation, confidence interval and p-value}
  \item{est.spillover}{Information of spillover effect estimation, containing point estimation, confidence interval and p-value}
}
\references{Brodersen et al. \emph{Inferring causal impact using Bayesian structural time-series models}. Annals of Applied Statistics, 2015}
\author{Zihao Zheng and Feiyu Yue}
\seealso{See also \code{?ad_cost}}
\examples{
  ## simulate data
  
  set.seed(1)

  y0 = 100 + arima.sim(model = list(ar = 0.3), n = 125)

  y.con = y0 + rnorm(125)
  y.con[101:125] = y.con[101:125] - 10 ## -10 as spillover effect

  y.exp = y0 + rnorm(125)
  y.exp[101:125] = y.exp[101:125] + 10 ## 10 as treatment effect

  pre.period = c(1:100)
  post.period = c(101:125)

  ## visualize

  plot(y.exp, col = "red", type = "l", ylab = "response",
     ylim = c(80, 120))

  lines(y.con, col = "blue")

  abline(v = 101, col = "grey", lty = 2, lwd = 2)

  legend("topleft", legend = c("exp", "con"), col = c("red", "blue"),
       cex = 1, lty = 1)

  ## try SPORTSCausal with ARIMA + AIC

  fit.aic = sportscausal(y.exp = y.exp, y.con = y.con, 
            pre.period = pre.period, post.period = post.period, is.plot = FALSE)

  fit.aic$est.treatment
  fit.aic$est.spillover

  ## you can also try model.select = "CV" or "lstm"
}

\keyword{Causal Inference}
\keyword{Spillover effect}
\keyword{Time series forecasting}
\keyword{Randomized Controlled Trial}
