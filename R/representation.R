# trim only for relationship, cutoff can apply in general?
#' Extract first kasumi representation
#'
#' Importance of relationships or WCC composition
#'
#' @export
extract_representation <- function(kasumi.results, cutoff = 0, trim = 1, ...) {
  repr.type <- ifelse(tibble::is_tibble(kasumi.results),
                      "composition", "relationship"
  )

  if (repr.type == "relationship") {
    sig <- extract_signature(kasumi.results,
                             type = "i",
                             intersect.targets = FALSE, trim = trim,
                             ...
    )
    sig[is.na(sig)] <- floor(min(sig %>% dplyr::select(-sample), na.rm = TRUE))
    sig <- sig %>% dplyr::mutate(dplyr::across(
      !sample,
      ~ ifelse(.x <= cutoff, 0, .x)
    ))

    keep <- which(sig %>% dplyr::select(-sample, -dplyr::contains("intra_")) %>%
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


#' Extract cluster based representation
#'
#' @export
extract_clusters <- function(kasumi.representation,
                             type = c("leiden", "kmeans"), sim = 0.8, res = 0.8,
                             k = 10, seed = 1) {
  clus.type <- match.arg(type)

  if (clus.type == "leiden") {
    clusters <- leiden(
      kasumi.representation %>%
        dplyr::select(-c(id, xcenter, ycenter)),
      sim, res, "cosine", seed
    )
  } else {
    clusters <- kmeans(
      kasumi.representation %>%
        dplyr::select(-c(id, xcenter, ycenter)),
      k, seed
    )
  }

  cluster.repr <- purrr::map(clusters, ~ .x == seq(length(unique(clusters)))) %>%
    purrr::reduce(rbind) %>%
    tibble::as_tibble(.name_repair = "unique") %>%
    dplyr::rename_with(~ stringr::str_remove_all(.x, "\\."))

  cbind(
    kasumi.representation %>%
      dplyr::select(id, xcenter, ycenter),
    cluster.repr
  )
}


#' Identify persistent clusters
#'
#' Representation is either raw kasumi.clusters or aggregated kasumi.clusters
#' The parameter is either a percentage if < 1 or a minimum number of samples
#' when > 1.
#'
#' @export
persistent_clusters <- function(kasumi.clusters, parameter = 0.1) {
  per.sample <- kasumi.clusters %>% dplyr::select(-id)

  if (length(unique(kasumi.clusters %>% dplyr::pull(id))) < nrow(kasumi.clusters)) {
    per.sample <- kasumi.clusters %>%
      dplyr::select(-c(xcenter, ycenter)) %>%
      dplyr::group_by(id) %>%
      dplyr::group_split(.keep = FALSE) %>%
      purrr::map_dfr(~ colSums(.x) / sum(.x))
  }

  par <- ifelse(parameter <= 1,
                parameter * nrow(per.sample),
                parameter
  )

  per.sample %>%
    dplyr::select(dplyr::where(~ ((sd(.) > 1e-3) | is.na(sd(.))) &
                                 (sum(. > 0) >= par))) %>%
    colnames()
}

#' @noRd
freq_repr <- function(labels) {
  colSums(labels) / sum(labels)
}

#' Represent each sample by its cluster composition
#' @param kasumi.clusters clusters calculated from the extracted representation
#'
#' @export
aggregate_clusters <- function(kasumi.clusters) {
  # originally persistent clusters were removed after distribution
  group.map <- kasumi.clusters %>%
    dplyr::select(-c(xcenter, ycenter)) %>%
    dplyr::group_by(id) %>%
    dplyr::group_split()

  group.ids <- group.map %>% purrr::map_chr(~ .x$id[1])

  group.map %>% purrr::map_dfr(~ .x %>% select(-id) %>% freq_repr()) %>%
    tibble::add_column(
      id = group.ids,
      .before = 1
    )
}
