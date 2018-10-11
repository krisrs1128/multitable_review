#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Helper functions for plotting multitable scores
##
## author: sankaran.kris@gmail.com
## date: 09/26/2017

library("ggplot2")
library("ggrepel")

#' Combine Loadings
#'
#' This combines loadings across data types. It's just a helper for the
#' WELL-China dataset
prepare_loadings <- function(loadings_list, types, K = 3) {
  df_list <- list()
  for (m in seq_along(loadings_list)) {
    colnames(loadings_list[[m]]) <- NULL
    df_list[[m]] <- data.frame(
      "variable" = rownames(loadings_list[[m]]),
      "seq_num" = NA,
      "type" = types[m],
      "Axis" = loadings_list[[m]][, 1:K]
    )
  }

  seq_ix <- which(types == "seq")
  if (length(seq_ix) > 0) {
    df_list[[seq_ix]]$seq_num <- df_list[[seq_ix]]$variable
  }

  do.call(rbind, df_list)
}

prepare_scores <- function(scores_list, types, K = 3) {
  df_list <- list()
  for (m in seq_along(scores_list)) {
    colnames(scores_list[[m]]) <- NULL
    df_list[[m]] <- data.frame(
      "number" = rownames(scores_list[[m]]),
      "type" = types[m],
      "Axis" = scores_list[[m]][, 1:K]
    )
  }

  do.call(rbind, df_list)
}

melt_scores <- function(scores) {
  scores %>%
    select(number, type, starts_with("Axis")) %>%
    gather(comp, value, starts_with("Axis")) %>%
    unite(comp_type, comp, type) %>%
    spread(comp_type, value)
}

reshape_posterior_score <- function(xi, bc) {
  xi %>%
    mutate(
      number = bc$number[row],
      col = paste0("axis_", col)
    ) %>%
    spread(col, value) %>%
    left_join(bc)
}

perc_label <- function(eigs, i) {
  perc <- 100 * eigs[i] / sum(eigs)
  sprintf("Axis %s [%s%%]", i, round(perc, 1))
}

taxa_order <- function(loadings) {
  hc <- loadings %>%
    select(starts_with("Axis")) %>%
    as.matrix() %>%
    dist() %>%
    hclust()

  hc$order
}

plot_topics <- function(loadings) {
  taxa_hc <- taxa_order(loadings)
  mloadings <- loadings %>%
    filter(type == "seq") %>%
    select(variable, starts_with("Axis"), family) %>%
    gather(topic, loading, -variable, -family)

  mloadings$variable <- factor(
    mloadings$variable,
    levels = loadings$seq_num[taxa_hc]
  )
  mloadings$topic <- gsub("Axis\\.", "Axis ", mloadings$topic)

  ggplot(mloadings) +
    geom_hline(yintercept = 0) +
    geom_point(
      aes(
        x = variable,
        col = family,
        y = loading
      )
    ) +
    facet_grid(topic ~ family, scale = "free_x", space = "free") +
    theme(
      panel.spacing.x = unit(0, "cm"),
      axis.text.x = element_blank(),
      strip.text.x = element_blank()
    )
}

plot_loadings <- function(loadings,
                          eigs,
                          size_breaks = c(-5, 5),
                          a = 1,
                          plot_dims = c(1, 2, 3)) {
  ggplot(loadings) +
    geom_hline(yintercept = 0, size = 0.5) +
    geom_vline(xintercept = 0, size = 0.5) +
    geom_point(
      data = loadings %>%
        filter(type == "seq"),
      aes_string(
        x = paste0("Axis.", plot_dims[1]),
        y = paste0("Axis.", plot_dims[2]),
        size = paste0("Axis.", plot_dims[3]),
        col = "family"
      ),
      alpha = a
    ) +
    geom_text_repel(
      data = loadings %>%
        filter(type == "body_comp"),
      aes_string(
        x = paste0("Axis.", plot_dims[1]),
        y = paste0("Axis.", plot_dims[2]),
        size = paste0("Axis.", plot_dims[3]),
        label = "variable"
      ),
      segment.size = 0.3,
      segment.alpha = 0.5,
      force = 0.05
    ) +
    labs(
      "x" = perc_label(eigs, plot_dims[1]),
      "y" = perc_label(eigs, plot_dims[2]),
      "size" = perc_label(eigs, plot_dims[3]),
      "col" = "Family"
    ) +
    scale_size_continuous(
      range = c(0.2, 2),
      breaks = size_breaks
    ) +
    guides(col=guide_legend(keyheight=0.1, default.unit="inch"),
           size=guide_legend(keyheight=0.1, default.unit="inch")) +
    coord_fixed(sqrt(eigs[plot_dims[2]] / eigs[plot_dims[1]]))
}

plot_scores <- function(scores, col_var, col_label, eigs, size_breaks = c(-3, 3)) {
  ggplot() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    geom_point(
      data = scores,
      aes_string(
        x = "Axis.1",
        y = "Axis.2",
        size = "Axis.3",
        col = col_var
      )
    ) +
    labs(
      "col" = col_label,
      "x" = perc_label(eigs, 1),
      "y" = perc_label(eigs, 2),
      "size" = perc_label(eigs, 3)
    ) +
    scale_size_continuous(range = c(0, 1.5), breaks = size_breaks) +
    coord_fixed(sqrt(eigs[2] / eigs[1]))
}

link_scores <- function(mscores, alpha = 0.1) {
  geom_segment(
    data = mscores,
    aes(
      x = Axis.1_body_comp, xend = Axis.1_seq,
      y = Axis.2_body_comp, yend = Axis.2_seq,
      size = (Axis.3_body_comp + Axis.3_seq) / 2
    ),
    alpha = alpha
  )
}

#' Useful for plotting scores for several features
plot_scores_wrapper <- function(xi, raw, processed, scv) {
  xi_df <- data.frame(
    "Axis" = xi,
    "number" = rownames(processed$bc)
  ) %>%
    left_join(raw$bc) %>%
    left_join(family_means(processed$mseqtab))
  list(
    plot_scores(xi_df, "age", "age", c(1, 1)) + scv,
    plot_scores(xi_df, "bmi", "BMI", c(1, 1)) + scv,
    plot_scores(xi_df, "Trunk_LM", "Trunk LM", c(1, 1)) + scv,
    plot_scores(xi_df, "rl_ratio", "Rum. / Lachn.", c(1, 1)) + scv
  )
}

melt_parameters <- function(theta_samples) {
  mtheta <- list()
  for (i in seq_along(theta_samples)) {
    if (length(dim(theta_samples[[i]])) > 2) {
      mtheta[[i]] <- theta_samples[[i]] %>%
        melt(
          varnames = c("iteration", "row", "col")
        )
    } else {
      mtheta[[i]] <- theta_samples[[i]] %>%
        melt(value.name = names(theta_samples)[i])
    }
  }
  names(mtheta) <- names(theta_samples)

  mtheta
}

seq_families <- function(mseqtab) {
  sf <- mseqtab %>%
    select(seq_num, family) %>%
    unique()
  sf$family <- factor(
    sf$family,
    names(sort(table(sf$family), decreasing = TRUE))
  )
  sf
}

mass_ordering <- function() {
  site_ordered <- c(
    "aoi", "fat_lean_ratio", "age", "height_dxa", "weight_dxa", "bmi",
    "android_fm", "android_lm", "gynoid_fm", "gynoid_lm", "l_trunk_fm",
    "l_trunk_lm", "r_trunk_fm", "r_trunk_lm", "trunk_fm", "trunk_lm",
    "l_total_fm", "l_total_lm", "r_total_fm", "r_total_lm", "total_fm",
    "total_lm", "l_leg_fm", "l_leg_lm", "r_leg_fm", "r_leg_lm", "legs_fm",
    "legs_lm", "l_arm_fm", "l_arm_lm", "r_arm_fm", "r_arm_lm", "arms_fm",
    "arms_lm"
  )
  c(
    site_ordered[!grepl("fm|lm", site_ordered)],
    site_ordered[grepl("fm", site_ordered)],
    site_ordered[grepl("lm", site_ordered)]
  )
}
