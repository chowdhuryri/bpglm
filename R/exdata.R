#' Bivariate count data set on 5567 respondents.
#'
#' A dataset containing the number of conditions and number of health care services utilizations during two years periods.
#' The data set comes from Health and Retirement Study (HRS) in US.
#' Full data set is freely availabible from HRS site after registration.
#' First two columns are two outcomes (y1, y2). Columns 3 to 6 are four covariates.
#' @format A data frame with 5567 rows and 10 variables:
#' \describe{
#'   \item{r10conde}{Number of conditions during a follow-up}
#'   \item{utiliza10}{Utilization of various health care services}
#'   \item{Gender}{Gender of the subject. Male=1; Female=0.}
#'   \item{Age}{Age of respondent (in years)}
#'   \item{Hispanic}{Race of the respondents. Hispanic=1; Others=0.}
#'   \item{Veteran}{Veteran status. Yes=1; No=0.}
#' }
#' @source \url{http://hrsonline.isr.umich.edu/}
"exdata"

