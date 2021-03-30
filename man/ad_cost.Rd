\name{ad_cost}
\alias{ad_cost}
\docType{data}
\title{Advertising cost: a real experimental data under spillover effect}
\description{This dataset comes from an A/Btest, which is to evaluate how a newly-proposed algorithm will affect the cost of advertising. Assuming that the bidding environment of an advertising market is stable in a short period of time, there will be no net increase or decrease of cost. When a treatment is applied, the mutual interference between the experiment group and the control group can not be ignored. For example, the difference in cost between the experiment group and the control group might not only come from the increase of experiment group, caused by treatment effect, but also from the potential decrease in control group. That is the typical situation for spillover causal inference to be implemented.}
\usage{data("ad_cost")}
\format{
  A data frame with 49 observations on the following 3 variables.
  \describe{
    \item{\code{y.exp}}{A numeric vector of responses in experiment group.}
    \item{\code{y.con}}{A numeric vector of responses in control group.}
    \item{\code{time}}{A numeric vector indicating time period before/after the treatment, \code{time = 1} represents post treatment period.}
  }
}
\details{This data has been linearly transformed for confidential issue.}
\examples{
### load data

data(ad_cost)

### define variables and visualize

y.exp = ad_cost$y.exp

y.con = ad_cost$y.con

plot(y.exp, col = "red", type = "l", 
     xlab = "time", ylab = "response")

lines(y.con, col = "blue")

### fit the model and return treatment/spillover effect

# notice that day-34 is the first day of treatment

fit = sportscausal(y.exp = y.exp, y.con = y.con,
                   pre.period = c(1:33), post.period = c(34:49), is.plot = FALSE)

fit$est.treatment

fit$est.spillover
}
\keyword{datasets}
