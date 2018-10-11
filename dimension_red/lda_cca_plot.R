#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Not generic enough to be included in plot.R, but don't want to copy code
## everywhere.
##
## author: sankaran.kris@gmail.com
## date: //2017

lda_cca_plots <- function(mdist, seq_fam, processed, opts) {
  scv <- scale_color_viridis(
    guide = guide_colorbar(barheight = 0.15, ticks = FALSE)
  )
  ggplot() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    geom_point(
      data = mdist$xi_s,
      aes(
        x = axis_1,
        y = axis_2,
        col = android_fm
      ),
      size = 0.1,
      alpha = 0.01
    ) +
    coord_equal() +
    scv +
    labs(x = "Axis 1", y = "Axis 2", col = "android_fm")

  ggsave(
    sprintf("%s/shared_scores_fm_posterior.png", opts$outdir),
    width = 5.63, height = 3.42
  )

  ## using the scores from just the bc table
  ggplot() +
    geom_hline(yintercept = 0, alpha = 0.5) +
    geom_vline(xintercept = 0, alpha = 0.5) +
    geom_point(
      data = mdist$xi_y,
      aes(
        x = axis_1,
        y = axis_2,
        col = android_fm
      ),
      size = 0.1,
      alpha = 0.1
    ) +
    coord_equal() +
    scv +
    labs(x = "Axis 1", y = "Axis 2", col = "Android FM")
  ggsave(
    sprintf("%s/unshared_scores_fm_posterior.png", opts$outdir),
    width = 9.24,
    height = 2.5
  )

  ## Now plotting loadings boxplots
  mdist$Wx$seq_num <- taxa_names(processed$ps)[mdist$Wx$row]
  mdist$Wx <- mdist$Wx %>%
    left_join(seq_fam)
  mdist$Wx$seq_num <- factor(
    mdist$Wx$seq_num,
    levels = names(sort(taxa_sums(processed$ps), decreasing = TRUE))
  )

  wx_summary <- mdist$Wx %>%
    group_by(col, seq_num) %>%
    dplyr::summarise(
      family = family[1],
      lower = quantile(value, 0.25),
      upper = quantile(value, 0.75),
      med = median(value)
    )

  ggplot(wx_summary) +
    geom_hline(yintercept = 0, alpha = 0.4) +
    geom_pointrange(
      aes(x = seq_num, ymin = lower, ymax = upper, y = med, col = family, fill = family),
      fatten = 1.2,
      size = 0.1
    ) +
    facet_grid(col ~ family, scale = "free_x", space = "free_x") +
    theme(
      panel.spacing.x = unit(0, "cm"),
      axis.text.x = element_blank(),
      strip.text.x = element_blank(),
      legend.position = "bottom"
    )
  ggsave(
    sprintf("%s/within_loadings_seq_boxplots.png", opts$outdir),
    width = 7.45,
    height = 3.37
  )

  mdist$Bx$seq_num <- taxa_names(processed$ps)[mdist$Bx$row]
  mdist$Bx <- mdist$Bx %>%
    left_join(seq_fam)
  mdist$Bx$seq_num <- factor(
    mdist$Bx$seq_num,
    levels = names(sort(taxa_sums(processed$ps), decreasing = TRUE))
  )

  bx_summary <- mdist$Bx %>%
    group_by(col, seq_num) %>%
    dplyr::summarise(
      family = family[1],
      lower = quantile(value, 0.25),
      upper = quantile(value, 0.75),
      med = median(value)
    )

  ggplot(bx_summary) +
    geom_hline(yintercept = 0, alpha = 0.4) +
    geom_pointrange(
      aes(x = seq_num, ymin = lower, ymax = upper, y = med, col = family, fill = family),
      fatten = 1.2,
      size = 0.1
    ) +
    facet_grid(col ~ family, scale = "free_x", space = "free_x") +
    theme(
      panel.spacing.x = unit(0, "cm"),
      axis.text.x = element_blank(),
      strip.text.x = element_blank(),
      legend.position = "bottom"
    )
  ggsave(
    sprintf("%s/between_loadings_seq_boxplots.png", opts$outdir),
    width = 7.45,
    height = 3.37
  )

  ## and finally loadings boxplot for body composition variables
  bc_names <- processed$ps %>%
    sample_data() %>%
    select(-id, -number, -gender, -batch, -operator) %>%
    colnames()
  mdist$Wy$variable <- factor(
    bc_names[mdist$Wy$row],
    levels = mass_ordering()
  )

  ggplot(mdist$Wy) +
    geom_hline(yintercept = 0, alpha = 0.4) +
    geom_boxplot(
      aes(
        x = variable,
        y = value
      ),
      outlier.size = 1
    ) +
    facet_grid(col ~ .) +
    theme(axis.text.x = element_text(angle = -90))

  ggsave(
    sprintf("%s/within_loadings_body_comp_boxplots.png", opts$outdir),
    width = 6.44,
    height = 3.95
  )

  mdist$By$variable <- factor(
    bc_names[mdist$By$row],
    levels = mass_ordering()
  )

  ggplot(mdist$By) +
    geom_hline(yintercept = 0, alpha = 0.4) +
    geom_boxplot(
      aes(
        x = variable,
        y = value
      ),
      outlier.size = 1
    ) +
    facet_grid(col ~ .) +
    theme(axis.text.x = element_text(angle = -90))

  ggsave(
    sprintf("%s/between_loadings_body_comp_boxplots.png", opts$outdir),
    width = 6.44,
    height = 3.95
  )
}
