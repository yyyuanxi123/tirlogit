#' Synthetic stated preference data for healthcare interventions
#'
#' @description A generated dataset `health` containing simulated choices and priority rankings
#' for medical interventions with different attribute, allowing for ties.
#'
#' The synthetic data includes three alternatives (A, B, and C), each representing a distinct medical
#' intervention. The profile of each intervention is a combination of seven attributes:
#' the type, the success probability (i.e. effectiveness), and the frequency of adverse effects,
#' of the intervention; severity of illness; patient responsibility for causing disease;
#' time span until symptoms appear; age at onset of illness.
#'
#' We simulate a scenario with 100 individuals, each receiving 5 choice tasks.
#' For each task, a priority ranking is generated for the interventions, where a lower value in rank indicates a higher priority.
#' For instance, rank 1 represents the highest priority, rank 2 a slightly lower priority, with larger values indicating progressively lower priorities.
#' Additionally, choice data is generated to indicate which intervention(s) is/are preferred to be reimbursed,
#' allowing for multiple best choices in cases where respondents have no strict preference.
#'
#' @name health
#' @docType data
#' @format A data frame in wide format with the following columns:
#' \describe{
#'   \item{id}{individual ID.}
#'   \item{taskid}{choice task ID.}
#'   \item{type.A, type.B, type.C}{type of the medical intervention for
#'   alternatives A, B, and C (\code{1} = preventive, \code{2} = curative).}
#'   \item{timeSpan.A, timeSpan.B, timeSpan.C}{time span until symptoms appear,
#'   (\code{1} = after 20 years, \code{2} = after 5 years, \code{3} = within a year). }
#'   \item{effectivenss.A, effectivenss.B, effectivenss.C}{the success probability (i.e. effectiveness) of the intervention
#'   (\code{1} = 33% effective, \code{2} = 66% effective, \code{3} = 100% effective).}
#'   \item{adverseEffects.A, adverseEffects.B, adverseEffects.C}{the frequency of adverse effects
#'   (\code{1} = often, \code{2} = rarely, \code{3} = never).}
#'   \item{severity.A, severity.B, severity.C}{severity of illness
#'   (\code{1} = not severe, \code{2} = severe, \code{3} = lethal).}
#'   \item{lifestyle.A, lifestyle.B, lifestyle.C}{responsibility of patient's lifestyle for causing disease,
#'   (\code{1} = fully, \code{2} = partly, \code{3} = not at all).}
#'   \item{ageGroup.A, ageGroup.B, ageGroup.C}{age at onset of illness
#'   (\code{1} = aged 80-90, \code{2} = aged 60-70, \code{3} = aged 40-50,
#'   \code{4} = aged 20-30, \code{5} = aged 0-10).}
#'   \item{rank.A, rank.B, rank.C}{reimbursement priority ranking for alternatives A, B, and C (lower values indicate higher priority).}
#'   \item{best.A, best.B, best.C}{an indicator variable representing whether intervention(s) are considered for reimbursement (\code{1} = yes, \code{2} = no).}
#' }
#'
#' @usage data("health", package = "groltir")
"health"


#' Synthetic stated preference data for transportation modes
#'
#' @description A generated dataset containing simulated choice data
#' for transportation modes across multiple choice tasks.
#' The dataset is structured in long format, where each row represents a specific transportation alternative
#' (e.g., train, air, bus, and car) within a choice task for an individual. Each alternative is characterized
#' by a set of attributes, including cost, in-vehicle and out-of-vehicle time, and frequency.
#'
#' The data includes generated responses across five different types of choice tasks types, with responses in specific variables:
#' - **Best choice (`best`)**: Indicate the most preferred alternative in each task, allowing only a single most preferred choice per task.
#' - **Best choice(s) (`multiBest`)**: Indicate the most preferred alternative in each task, allowing for multiple best choices.
#' - **Ranking not allowing for ties (`rankUnique`)**: Provide a preference ranking for each alternative within a task, without allowing for ties.
#' - **Ranking allowing for ties (`rank`)**: Provide a preference ranking that allows ties, where lower values indicate higher preference. Tied alternatives share the same rank.
#' - **Best-worst choice(s) (`BestWorst`)**: Indicate the most preferred and the least preferred alternative in each task.
#'
#' @name mode_long
#' @docType data
#' @format A data frame in long format with the following columns:
#' \describe{
#'   \item{case}{individual ID.}
#'   \item{taskid}{choice task ID.}
#'   \item{alt}{name of the alternative (i.e. transportation mode, one of train, car, bus and air).}
#'   \item{dist}{choice set specific variable, the distance of the trip.}
#'   \item{cost}{alternative specific variable, monetary cost associated with the transportation mode.}
#'   \item{ivt}{alternative specific variable, in-vehicle travel time.}
#'   \item{ovt}{alternative specific variable, out-of-vehicle travel time.}
#'   \item{freq}{alternative specific variable, frequency of availability for the transportation mode.}
#'   \item{income}{choice set specific variable, household income.}
#'   \item{urban}{choice set specific variable, a dummy variable indicating whether the trip involves a large city as either
#'   the origin or the destination (\code{1} = yes, \code{0} = no).}
#'   \item{best}{a dummy variable indicating whether the alternative is selected as the best alternative; only one alternative can be selected as the best
#'   in each choice task (\code{1} = most preferred, \code{0} = otherwise).}
#'   \item{multiBest}{a dummy variable indicating whether the alternative is selected as the best alternative; multiple best choices can be given (\code{1} = most preferred, \code{0} = otherwise).}
#'   \item{rankUnique}{preference ranking of the alternative; each alternative is assigned a unique rank; rank 1 for the most preferred alternative.}
#'   \item{rank}{preference ranking of the alternative; ties are allowed; rank 1 for the most preferred alternative(s).}
#'   \item{BestWorst}{an indicator variable representing whether the alternative is selected as the most preferred or least preferred alternative; multiple best and worst choices are allowed (\code{1} = most preferred, \code{-1} = least preferred, \code{0} = otherwise).}
#' }
#' @source This dataset is generated based on the `ModeCanada` dataset from the \code{R} package `mlogit`.
#' @references
#' Croissant, Y. (2020). `mlogit`: Multinomial Logit Models. R package version 1.1-1.
#' Available from: \url{https://CRAN.R-project.org/package=mlogit}.
#' @usage data(mode_long, package = "groltir")
"mode_long"

#' Synthetic stated preference data for transportation modes
#'
#' @description A generated dataset containing simulated choice data
#' for transportation modes across multiple choice tasks, in a wide format. This dataset is a wide-format
#' version of \link{mode_long}. Each row corresponds to a choice task.
#' It includes choice set specific variables (`dist`, `income`, and `urban`),
#' alternative specific attributes for each mode (`cost`, `freq`, `ivt`, and `ovt`),
#' and simulated preference data across various scenarios in separate columns.
#' For detailed descriptions of the dataset, see \link{mode_long}.
#' Additionally, this dataset also includes character string representations of each preference response.
#'
#' @name mode_wide
#' @docType data
#' @format A data frame in wide format.
#' @usage data(mode_wide, package = "groltir")
"mode_wide"
