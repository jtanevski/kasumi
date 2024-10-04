# Kasumi runner
# Copyleft (ɔ) 2024 Jovan Tanevski [jovan.tanevski@uni-heidelberg.de]

#' @importFrom rlang !! :=
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Kasumi is able to run computationally intensive functions
  in parallel. Please consider specifying a future::plan(). For example by running
  future::plan(future::multisession) before calling Kasumi functions.")
}


#' Train Kasumi models
#'
#' Train local multi-view models by sliding a window across the sample as captured by the
#' view composition.
#'
#' @param views view composition.
#' @param positions a \code{data.frame}, \code{tibble} or a \code{matrix}
#'     with named coordinates in columns and rows for each spatial unit ordered
#'     as in the intraview.
#' @param window size of the window.
#' @param overlap overlap of consecutive windows (percentage).
#' @param sample.id id of the sample.
#' @param results.db path to the database file to store the results.
#' @param minu minimum number of spatial units in the window.
#' @param ... all other parameters are passed to \code{\link{run_misty}()}.
#'
#' @return Path to the result folder(s) that can be passed to
#'     \code{\link{collect_results}()}.
#'
#' @examples
#' # Create a view composition of an intraview and a paraview with radius 10 then
#' # run Kasumi for a single sample.
#'
#' library(dplyr)
#'
#' # get the expression data
#' data("synthetic", package = "mistyR")
#' expr <- synthetic[[1]] %>% select(-c(row, col, type))
#' # get the coordinates for each cell
#' pos <- synthetic[[1]] %>% select(row, col)
#'
#' # compose
#' kasumi.views <- create_initial_view(expr) %>% add_paraview(pos, l = 10)
#'
#' # run with a window of size 100
#' run_kasumi(kasumi.views, pos, window = 100)
#'
#' @export
run_kasumi <- function(views, positions, window, overlap = 50,
                       sample.id = "sample",
                       results.db = paste0(sample.id, ".sqm"), minu = 50,
                       ...) {
  db.file <- R.utils::getAbsolutePath(results.db)
  db.lock <- paste0(db.file, ".lock")

  if (!file.exists(db.file)) {
    current.lock <- filelock::lock(db.lock)
    mistyR:::create_sqm(db.file)
    filelock::unlock(current.lock)
    rm(current.lock)
  }

  current.lock <- filelock::lock(db.lock)
  sqm <- DBI::dbConnect(RSQLite::SQLite(), db.file)
  if (!DBI::dbExistsTable(sqm, "wcc")) {
    DBI::dbCreateTable(
      sqm, "wcc",
      c(sample = "TEXT", Target = "TEXT", value = "REAL")
    )
  }
  DBI::dbDisconnect(sqm)
  filelock::unlock(current.lock)

  x <- tibble::tibble(
    xl = seq(
      min(positions[, 1]),
      max(positions[, 1]),
      window - window * overlap / 100
    ),
    xu = xl + window
  ) %>%
    dplyr::filter(xl < max(positions[, 1])) %>%
    dplyr::mutate(xu = pmin(xu, max(positions[, 1]))) %>%
    round(2)

  y <- tibble::tibble(
    yl = seq(
      min(positions[, 2]),
      max(positions[, 2]),
      window - window * overlap / 100
    ),
    yu = yl + window
  ) %>%
    dplyr::filter(yl < max(positions[, 2])) %>%
    dplyr::mutate(yu = pmin(yu, max(positions[, 2]))) %>%
    round(2)

  tiles <- tidyr::expand_grid(x, y)

  # make nested plan here with sequential at the second level
  # retrieve current plan by simply running plan() without parameters
  old.plan <- future::plan()
  future::plan(list(old.plan, future::sequential))
  message("\nSliding")
  tiles %>% furrr::future_pwalk(\(xl, xu, yl, yu){
    selected.rows <- which(
      positions[, 1] >= xl & positions[, 1] <= xu &
        positions[, 2] >= yl & positions[, 2] <= yu
    )

    if (length(selected.rows) >= minu) {
      filtered.views <- views %>%
        filter_views(selected.rows)

      suppressMessages(mistyR::run_misty(
        filtered.views,
        paste0(sample.id, "/", xl, "_", yl, "_", xu, "_", yu),
        results.db,
        ...
      ))

      wcc(
        filtered.views[["intraview"]],
        paste0(sample.id, "/", xl, "_", yl, "_", xu, "_", yu),
        results.db
      )
    }
  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))

  future::plan(old.plan)

  if (file.exists(db.lock)) file.remove(db.lock)
  return(results.db)
}


#' Calculate window composition and store
#' @noRd
wcc <- function(intra.view, sample.id, results.db) {
  composition <- (intra.view %>% colSums()) / sum(intra.view)

  to.write <- t(composition) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      tidyr::everything(),
      names_to = "Target",
      values_to = "value"
    ) %>%
    tibble::add_column(sample = sample.id, .before = 1)

  db.file <- R.utils::getAbsolutePath(results.db)
  db.lock <- paste0(db.file, ".lock")

  current.lock <- filelock::lock(db.lock)
  sqm <- DBI::dbConnect(RSQLite::SQLite(), db.file)
  DBI::dbAppendTable(sqm, "wcc", to.write)
  DBI::dbDisconnect(sqm)
  filelock::unlock(current.lock)
}
