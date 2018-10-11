#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Sparse partial least squares applied to the WELL microbiome data.
##
## author: sankaran.kris@gmail.com
## date: 10/12/2017

###############################################################################
## Libraries and setup
###############################################################################
library("phyloseq")
library("spls")
library("reshape2")
library("tidyverse")
source("../dimension_red/prep_tables.R")
source("../dimension_red/plot.R")
set.seed(20171101)

## cleaner ggplot theme
scale_colour_discrete <- function(...)
  scale_color_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6', "#464646"),
    na.value = "black"
  )
scale_fill_discrete <- function(...)
  scale_fill_manual(
    values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6', "#464646"),
    na.value = "black"
  )

theme_set(theme_bw())
theme_update(
  panel.background = element_rect(fill = "#F8F8F8"),
  panel.border = element_rect(size = 0.5),
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6),
  axis.text = element_text(size = 6),
  axis.title = element_text(size = 8),
  strip.background = element_blank(),
  strip.text = element_text(size = 8),
  legend.key = element_blank()
)

###############################################################################
## read and prepare the data
###############################################################################
raw <- read_data()
processed <- process_data(
  raw$seqtab,
  raw$bc,
  raw$bc_full,
  raw$taxa,
  raw$tree
)

y_df <- sample_data(processed$ps)
keep_ix <- apply(y_df, 1, function(x) all(!is.na(x)))
rm_cols <- which(colnames(y_df) %in% c("batch", "operator", "gender", "number", "id"))
y <- scale(y_df[keep_ix, -rm_cols])
x <- get_taxa(processed$ps)
cv_eval <- cv.spls(x, y, K = 4:8, eta = seq(0, 0.6, 0.05), scale.x = FALSE, fold = 5)
cv_eval

train_ix <- sample(1:nrow(x), 80)
fit <- spls(x[train_ix, ], y[train_ix, ], cv_eval$K.opt, cv_eval$eta.opt)
y_hat <- x %*% fit$betahat
plot(y[train_ix, 24], y_hat[train_ix, 24])
points(y[-train_ix, 24], y_hat[-train_ix, 24], col = "blue")
abline(a = 0, b = 1, col = "red")

## refit on full data
fit <- spls(x, y, cv_eval$K.opt, cv_eval$eta.opt)

###############################################################################
## plot fitted coefficients
###############################################################################
seq_fam <- seq_families(processed$mseqtab)
mass_type_ordered <- mass_ordering()

dend <- as.dendrogram(hclust(dist(coef(fit))))
hc <- as.hclust(reorder(dend, coef(fit), mean))
species_order <- colnames(x)[hc$order]

mbeta <- coef(fit) %>%
  melt(
    varnames = c("seq_num", "feature"),
    value.name = "coef"
  ) %>%
  left_join(seq_fam) %>%
  mutate(
    feature = factor(feature, mass_type_ordered),
    seq_num = factor(seq_num, species_order)
  )

ggplot(mbeta) +
  geom_tile(
    aes(x = seq_num, y = feature, fill = coef)
  ) +
  geom_rect(
    aes(col = family),
    fill = "transparent", size = 2,
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  scale_fill_gradient2(
    "Coef. ",
    guide = guide_colorbar(ticks = FALSE, barheight = 0.9),
    low = "#40004b",
    high = "#00441b"
  ) +
  guides(
    color = guide_legend(nrow = 4)
  ) +
  scale_x_discrete("Genus", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  facet_grid(. ~ family, scale = "free", space = "free") +
  theme(
    axis.text.x = element_blank(),
    panel.spacing = unit(0, "cm"),
    strip.text.x = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  "../chapter/figure/spls/coef_heatmap.png",
  width = 5.9,
  height = 4.1
)

large_species <- mbeta %>%
  filter(
    feature %in% c("android_fm", "gynoid_fm")
  ) %>%
  group_by(seq_num) %>%
  dplyr::mutate(norm = sqrt(sum(coef ^ 2))) %>%
  filter(norm > 0.05) %>%
  arrange(desc(norm))

mlarge_species <- melt(
  data.frame(
    "Number" = rownames(x),
    x[, unique(as.character(large_species$seq_num))],
    y[, c("android_fm", "gynoid_fm")]
  ),
  id.vars = c("Number", "android_fm", "gynoid_fm"),
  variable.name = "seq_num"
) %>%
  left_join(seq_fam)

mlarge_species$seq_num <- factor(
  mlarge_species$seq_num,
  large_species %>%
    filter(feature == "gynoid_fm") %>%
    arrange(desc(coef)) %>%
    .[["seq_num"]]
)

ggplot(mlarge_species) +
  geom_hline(yintercept = 0, size = 0.1, alpha = 0.8) +
  stat_smooth(
    aes(x = value, y = gynoid_fm, col = family),
    alpha = 0.8, size = 0.5, method = "lm", fill = "#dfdfdf"
  ) +
  geom_vline(xintercept = 0, size = 0.1, alpha = 0.8) +
  geom_point(
    aes(x = value, y = gynoid_fm, col = family),
    size = 0.7, alpha = 0.8
  ) +
  facet_wrap(~seq_num, ncol = 6) +
  theme(legend.position = "bottom")

ggsave(
  "../chapter/figure/spls/gynoid_fm_species.png",
  width = 7.4,
  height = 6.3
)

## same plot for total fm
mlarge_species$seq_num <- factor(
  mlarge_species$seq_num,
  large_species %>%
    filter(feature == "android_fm") %>%
    arrange(desc(coef)) %>%
    .[["seq_num"]]
)

ggplot(mlarge_species) +
  geom_hline(yintercept = 0, size = 0.1, alpha = 0.8) +
  geom_vline(xintercept = 0, size = 0.1, alpha = 0.8) +
  stat_smooth(
    aes(x = value, y = android_fm, col = family),
    alpha = 0.8, size = 0.5, method = "lm", fill = "#dfdfdf"
  ) +
  geom_point(
    aes(x = value, y = android_fm, col = family),
    size = 0.7, alpha = 0.8
  ) +
  facet_wrap(~seq_num, ncol = 6) +
  theme(legend.position = "bottom")

ggsave(
  "../chapter/figure/spls/android_fm_species.png",
  width = 7.4,
  height = 6.3
)
