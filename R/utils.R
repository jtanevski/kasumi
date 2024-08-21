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


# trim only for relationship, cutoff can apply in general?
#' Extract first kasumi representation
#'
#' Importance of relationships or KWC composition
#'
#' @export
extract_representation <- function(kasumi.results, cutoff = 0, trim = 1) {
  repr.type <- ifelse(tibble::is_tibble(kasumi.results),
    "composition", "relationship"
  )

  if (repr.type == "relationship") {
    sig <- extract_signature(kasumi.results,
      type = "i",
      intersect.targets = FALSE, trim = trim
    )
    sig[is.na(sig)] <- floor(min(sig %>% select(-sample), na.rm = TRUE))
    sig <- sig %>% dplyr::mutate(dplyr::across(
      !sample,
      ~ ifelse(.x <= cutoff, 0, .x)
    ))

    keep <- which(sig %>% dplyr::select(-sample, -contains("intra_")) %>%
      rowSums() != 0)

    samps <- sig %>%
      dplyr::slice(keep) %>%
      dplyr::select(sample) %>%
      dplyr::mutate(
        id = stringr::str_extract(sample, ".*/") %>% stringr::str_remove("/"),
        box = stringr::str_extract(sample, "/[0-9].*$") %>%
          stringr::str_remove("/")
      ) %>%
      dplyr::rowwise(id) %>%
      dplyr::summarize(rebox = box %>%
        stringr::str_split("_", simplify = T) %>%
        as.numeric() %>% list(), .groups = "drop") %>%
      dplyr::rowwise(id) %>%
      dplyr::summarize(
        xcenter = (rebox[3] + rebox[1]) / 2,
        ycenter = (rebox[4] + rebox[2]) / 2, .groups = "drop"
      )

    # the filtering here also matters
    clean <- sig %>%
      dplyr::select(-sample, -dplyr::contains("_.novar")) %>%
      dplyr::slice(keep) %>%
      dplyr::select(dplyr::where(~ sum(.) != 0))
  } else {
    wide <- kasumi.results %>% tidyr::pivot_wider(
      id_cols = "sample",
      names_from = "Target", values_from = "value"
    )

    samps <- wide %>%
      dplyr::select(sample) %>%
      dplyr::mutate(
        id = stringr::str_extract(sample, ".*/") %>% stringr::str_remove("/"),
        box = stringr::str_extract(sample, "/[0-9].*$") %>%
          stringr::str_remove("/")
      ) %>%
      dplyr::rowwise(id) %>%
      dplyr::summarize(rebox = box %>%
        stringr::str_split("_", simplify = T) %>%
        as.numeric() %>% list(), .groups = "drop") %>%
      dplyr::rowwise(id) %>%
      dplyr::summarize(
        xcenter = (rebox[3] + rebox[1]) / 2,
        ycenter = (rebox[4] + rebox[2]) / 2, .groups = "drop"
      )

    clean <- wide %>%
      dplyr::select(-sample) %>%
      dplyr::select(dplyr::where(~ sum(.) != 0))
  }

  return(cbind(samps, clean))
}



# Extract cluster based representation
#' @export
extract_clusters <- function(kasumi.representation,
                             type = c("leiden", "kmeans"), sim = 0.8, res = 0.8,
                             k = 10, persistence = 0.1, seed = 1) {
  clus.type <- match.arg(type)

  if (clus.type == "leiden") {
    clusters <- leiden(
      kasumi.representation %>%
        dplyr::select(-c(id, xcenter, ycenter)),
      sim, res, seed
    )
  } else {
    clusters <- kmeans(
      kasumi.representation %>%
        dplyr::select(-c(id, xcenter, ycenter)),
      k, seed
    )
  }

  cluster.repr <- purrr::map(clusters, ~ .x == seq(length(unique(clusters)))) %>%
    reduce(rbind) %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    dplyr::rename_with(~ stringr::str_remove_all(.x, "\\."))


  per <- ifelse(persistence <= 1,
    persistence * length(unique(kasumi.representation %>% pull(id))),
    persistence
  )

  persistent.clusters <- kasumi.representation %>%
    select(id) %>%
    cbind(cluster.repr) %>%
    dplyr::group_by(id) %>%
    dplyr::group_split(.keep = FALSE) %>%
    purrr::map_dfr(~ colSums(.x) / sum(.x)) %>%
    dplyr::select(dplyr::where(~ ((sd(.) > 1e-3) | is.na(sd(.))) & (sum(. > 0) >= per))) %>%
    colnames()

  cbind(kasumi.representation %>%
    dplyr::select(id, xcenter, ycenter),
    cluster.repr %>% dplyr::select(dplyr::all_of(persistent.clusters)))
}


# Represent each sample by its cluster composition
#' @export
aggregate_clusters <- function(kasumi.clusters) {
  #originally persistent clusters were removed after distribution
  kasumi.clusters %>% dplyr::select(-c(xcenter, ycenter)) %>%
    dplyr::group_by(id) %>% group_split() %>%
    purrr::map_dfr(~ colSums(.x %>% dplyr::select(-id)) / sum(.x %>% dplyr::select(-id))) %>%
    tibble::add_column(id = kasumi.clusters %>% pull(id) %>% unique(), .before = 1)
}

# helper functions for clustering

#' @noRd
leiden <- function(representation, minsim = 0.8, resolution = 0.8,
                   measure = "cosine", seed = 1) {
  sim <- proxy::simil(representation, measure)

  sim[sim < minsim] <- 0
  withr::with_seed(
    seed,
    groups <-
      igraph::graph.adjacency(sim %>% as.matrix(),
        mode = "undirected", weighted = TRUE
      ) %>%
      igraph::cluster_leiden(
        resolution_parameter = resolution,
        n_iterations = -1
      )
  )

  return(groups$membership)
}

#' @noRd
kmeans <- function(representation, k = 10, seed = 1) {
  withr::with_seed(
    seed,
    clust <- ClusterR::KMeans_rcpp(representation, k)
  )
  return(clust$clusters)
}
