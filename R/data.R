#' Multi-omics dataset of sleep deprivation in mouse
#'
#' @description A dataset containing transcript and metabolite values analysed
#' in Albanese et al. 2023, subset of the multi-omics data set published in
#' Jan, M., Gobet, N., Diessler, S. et al. A multi-omics digital research object
#' for the genetics of sleep regulation. Sci Data 6, 258 (2019).
#'
#' `multi_omics_sd_small` is a smaller version, limited to the transcript Cirbp
#' and the transcripts and metabolites belonging to its neighborhood as
#' described in Albanese et al. 2023
#'
#' `multi_omics_sd_micro` is a minimal version with Cirbp and a selection of its
#'  neighborhood.
#'
#' @format ## `multi_omics_sd` 
#' A data frame with 30 rows and 238 variables (162
#'   transcripts and 76 metabolites):
#' \describe{
#'   \item{Plin4 to Tfrc}{log2 CPM values of 162 transcripts in mouse cortex
#'    under sleep deprivation (-4.52--10.46)}
#'   \item{Ala to SM C24:1}{abundance values of 76 metabolites (0.02--1112.67)}
#' }
#' @source Jan, M., Gobet, N., Diessler, S. et al. A multi-omics digital
#'   research object for the genetics of sleep regulation. Sci Data 6, 258
#'   (2019) <https://doi.org/10.1038/s41597-019-0171-x>
#' @source Figshare folder of the original manuscript:
#'   <https://figshare.com/articles/dataset/Input_data_for_systems_genetics_of_sleep_regulation/7797434>
"multi_omics_sd"

#' @rdname multi_omics_sd
#' @format ## `multi_omics_sd_small`
#' A data frame with 30 rows and 21 variables (16 transcripts and 5 metabolites)
#' \describe{
#'   \item{Cirbp to Stip1}{log2 CPM values of 16 transcripts in mouse cortex
#'    under sleep deprivation (4.24--9.31)}
#'   \item{Phe to PC ae C32:2}{Abundance values of 5 metabolites (0.17--145.33)}
#' }
"multi_omics_sd_small"

#' @rdname multi_omics_sd
#' @format ## `multi_omics_sd_micro`
#' A data frame with 30 rows and 6 variables (4 transcripts and 2 metabolites)
#' \describe{
#'   \item{Cirbp to Dnajb11}{log2 CPM values of 4 transcripts in mouse cortex
#'    under sleep deprivation (4.78--9.31)}
#'   \item{Trp to PC aa C36:3}{Abundance values of 2 metabolites (58.80--145.33)}
#' }
"multi_omics_sd_micro"
