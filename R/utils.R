#' Collect KWC data
#'
#' Collect KWC composition data per window produced by \code{\link{run_kasumi}()}.
#'
#' @param db.file path to database file with raw results from
#'     \code{\link{run_kasumi}()}.
#' @param sample.pattern a regex pattern to match sample names.
#'
#' @return Long format \code{tibble} with the relative abundance of each
#'        \var{Target} in each window and \var{sample}.
#'
#' @export
collect_kwc <- function(db.file, sample.pattern = ".") {
  sqm <- DBI::dbConnect(RSQLite::SQLite(), db.file)

  message("\nCollecting KWC data")

  try.filter <- DBI::dbReadTable(sqm, "kwc") %>%
    tibble::as_tibble() %>%
    dplyr::filter(stringr::str_detect(sample, sample.pattern))

  if (nrow(try.filter) == 0) {
    warning("sample.pattern doesn't match anything in the results database")
    return(NULL)
  }

  DBI::dbDisconnect(sqm)

  return(try.filter)
}

cluster_kasumi <- function(kasumi.results) {}

cluster_kwc <- function(kwc.results) {}

# helper functions for leiden and kmeans
