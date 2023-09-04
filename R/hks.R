#' Replication data for Hultman, Kathman, and Shannon (2013) United Nations Peacekeeping and Civilian Protection in Civil War
#'
#' A reduced replication data set from \insertCite{hultman2013united;textual}{evi} United Nations Peacekeeping and Civilian Protection in Civil War. Used to reproduce the the results from \insertCite{randahl2023evzinb;textual}{evi}. Note, to reproduce any other results from \insertCite{hultman2013united;textual}{evi} please download the original replication dataset using the link under source.
#'
#' @format 
#' A tibble with 3746 rows and 12 columns:
#' \describe{
#'   \item{conflict_id}{The Uppsala Conflict Data Programme conflict ID for the conflict}
#'   \item{osvAll}{The number of observed fatalities from one-sided violence against civilian in the specified conflict-month}
#'   \item{troopLag}{The number of UN military troops in thousands of troops (lagged)}
#'   \item{policeLag}{The number of UN police in thousands of troops (lagged)}
#'   \item{militaryobserversLag}{The number of UN military troops in thousands of troops (lagged)}
#'   \item{epduration}{The number of months the current conflict-episode has been ongoing}
#'   \item{lntpop}{The natural logarithm of the population of the country in which the conflict takes place}
#'   \item{lnbrv_AllLag}{The natural logarithm of the total number of battle related deaths in the conflict in the previous month}
#'   \item{osvAllLagDum}{A dummy variable taking the value 1 if any one-sided violence against civilians took place in the previous conflict month}
#'   \item{incomp}{A dummy variable taking the value 1 if the conflict is about government and 0 otherwise}
#'   \item{lntroopLag}{The log1p logarithm of troopLag}
#'   \item{lnepdur}{The log1p logarithm of the episode duration}
#' }
#' @source https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/6EBCGA
#' @references 
#' \insertAllCited{}
"hks"