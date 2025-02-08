#' Collect WCC data
#'
#' Collect Window Composition Clustering result per window produced by \code{\link{run_kasumi}()}.
#'
#' @param db.file path to database file with raw results from
#'     \code{\link{run_kasumi}()}.
#' @param sample.pattern a regex pattern to match sample names.
#'
#' @return Long format \code{tibble} with the relative abundance of each
#'        \var{Target} in each window and \var{sample}.
#'
#' @export
collect_wcc <- function(db.file, sample.pattern = ".") {
  sqm <- DBI::dbConnect(RSQLite::SQLite(), db.file)

  message("\nCollecting WCC data")

  try.filter <- DBI::dbReadTable(sqm, "wcc") %>%
    tibble::as_tibble() %>%
    dplyr::filter(stringr::str_detect(sample, sample.pattern))

  if (nrow(try.filter) == 0) {
    warning("sample.pattern doesn't match anything in the results database")
    return(NULL)
  }

  DBI::dbDisconnect(sqm)

  return(try.filter)
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

#' Use Kasumi representation for downstream classification
#'
#' @export
downstream_classify <- function(kasumi.clusters, target, dir = "auto", seed = 1) {
  representation <- kasumi.clusters %>%
    dplyr::select(-id) %>%
    tibble::add_column(target = target)

  withr::with_seed(
    seed,
    suppressWarnings(
      model <- caret::train(target ~ ., representation,
        method = "glm", metric = "ROC",
        trControl = caret::trainControl(
          method = "cv", number = 10,
          classProbs = TRUE,
          summaryFunction = caret::twoClassSummary,
          savePredictions = TRUE
        )
      )
    )
  )

  pROC::roc(model$pred$obs, model$pred[, 3], direction = dir, quiet = TRUE)
}


#' Calculate signed Model Reliance
#'
#' @export
sMR <- function(kasumi.clusters, target, seed = 1) {
  representation <- kasumi.clusters %>%
    dplyr::select(-id) %>%
    tibble::add_column(target = as.factor(target))

  model <- stats::glm(target ~ ., representation, family = "binomial")

  eorig <- downstream_classify(kasumi.clusters, target, dir = "auto", seed)
  dir <- eorig$direction
  dir.sign <- ifelse(dir == ">", -1, 1)
  cat(paste0("AUC: ", eorig$auc))

  withr::with_seed(
    seed,
    splitr <- stats::runif(nrow(representation)) %>% rank()
  )

  nas <- names(which(is.na(stats::coef(model)[-1])))

  eswitch <- representation %>%
    dplyr::select(-target, -nas) %>%
    colnames() %>%
    purrr::map_dbl(\(cname){
      downstream_classify(
        kasumi.clusters %>% dplyr::select(-nas) %>%
          dplyr::mutate(!!cname := representation[splitr, cname] %>% unlist()),
        target, dir, seed
      )$auc
    })


  mr <- dir.sign * sign(stats::coef(model, complete = FALSE)[-1]) * (1 - eswitch) / (1 - eorig$auc)

  toreturn <- tibble::tibble(Cluster = as.factor(names(mr)), sMR = mr)

  print(ggplot2::ggplot(
    toreturn %>% dplyr::mutate(Cluster = forcats::fct_reorder(Cluster, sMR)),
    ggplot2::aes(x = Cluster, y = sMR)
  ) +
    ggplot2::geom_segment(ggplot2::aes(x = Cluster, xend = Cluster, y = 0, yend = sMR)) +
    ggplot2::geom_point(ggplot2::aes(x = Cluster, y = sMR, color = sMR)) +
    ggplot2::scale_color_steps2(low = "darkgreen", mid = "white", high = "blue3") +
    ggplot2::geom_hline(yintercept = 0, color = "gray50") +
    ggplot2::geom_hline(yintercept = 1, color = "gray70", linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -1, color = "gray70", linetype = "dashed") +
    ggplot2::geom_label(label = paste("\u2190", eorig$levels[2]), x = length(mr), y = -1) +
    ggplot2::geom_label(label = paste(eorig$levels[1], "\u2192"), x = 1, y = 1) +
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ))

  return(toreturn)
}

# from unaggregated kasumi.clusters
#' Collect Kasumi cluster description
#'
#' @export
collect_kasumi_cluster <- function(kasumi.clusters, cluster, db.file) {
  sm.repr.all <- kasumi.clusters

  left <- sm.repr.all %>%
    dplyr::filter(dplyr::if_any(!!cluster)) %>%
    dplyr::select(id, xcenter, ycenter)

  dbcon <- DBI::dbConnect(RSQLite::SQLite(), db.file)
  samples <- DBI::dbGetQuery(dbcon, "SELECT DISTINCT sample FROM contributions") %>%
    unlist()
  matching <- grep(paste0("(", paste0(unique(left$id), collapse = "|"), ")"),
    samples,
    value = TRUE
  )
  DBI::dbDisconnect(dbcon)

  right <- tibble::tibble(sample = matching) %>%
    dplyr::mutate(
      id = stringr::str_extract(sample, ".*/") %>% stringr::str_remove("/"),
      box = stringr::str_extract(sample, "/[0-9].*$") %>% stringr::str_remove("/")
    ) %>%
    dplyr::rowwise(id) %>%
    dplyr::summarize(sample = sample, rebox = box %>%
      stringr::str_split("_", simplify = T) %>%
      as.numeric() %>% list(), .groups = "drop") %>%
    dplyr::rowwise(id) %>%
    dplyr::summarize(
      sample = sample,
      xcenter = (rebox[3] + rebox[1]) / 2,
      ycenter = (rebox[4] + rebox[2]) / 2, .groups = "drop"
    )

  pattern <- paste0("(", paste0(left %>%
    dplyr::left_join(right, by = c("id", "xcenter", "ycenter")) %>%
    dplyr::pull(sample), collapse = "|"), ")")

  collect_results(db.file, pattern)
}
