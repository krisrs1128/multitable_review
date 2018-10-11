#! /usr/bin/env Rscript

## File description -------------------------------------------------------------
##
## Preprocessing steps to apply before applying dimensionality reduction
## methods.
##
## author: sankaran.kris@gmail.com
## date: 09/27/2017

library("phyloseq")
library("DESeq2")
library("tidyverse")
library("forcats")
library("reshape2")

###############################################################################
## Prepare opts and read data
###############################################################################
merge_process_opts <- function(opts = list()) {
  default_opts <- list(
    gender = "Female", ## has more data
    filter_k = 0.07,
    filter_a = 5,
    scale_sample_data = TRUE,
    transform_fun = identity,
    vst = TRUE,
    remove_na = FALSE,
    outdir = "../data"
  )
  modifyList(default_opts, opts)
}

#' Read raw data
#'
#' Reads and applies basic preprocessing steps to the raw data.
#'
#' @param data_dir [string] The path to the raw data, relative to the current
#'   working directory.
#' @return results [list] A list containing the following essentially raw data
#'   sets,
#'    (1) taxa: The table giving taxonomic information for each sequence.
#'    (2) seqtab: The raw counts table of species across samples.
#'    (3) bc: The body composition data on the subset of people in the
#'         microbiome study.
#'    (4) bc_full: The body composition data across all participants in the
#'         study.
read_data <- function(simulate=TRUE, data_dir = "../data/") {
  if (simulate) {
    return (get(load(file.path(data_dir, "sim.rda"))))
  }

  seqtab <- readRDS(file.path(data_dir, "seqtab.rds"))
  batch <- read_csv(file.path(data_dir, "222_Participant_Batch_Numbers.csv"))
  bc <- readRDS(file.path(data_dir, "sample_data_bc.rds")) %>%
    mutate(id = as.character(id)) %>%
    left_join(batch)
  colnames(bc) <- tolower(colnames(bc))

  colnames(seqtab) <- paste0("genus_", seq_len(ntaxa(seqtab)))
  bc_full <- read_csv(file.path(data_dir, "WELL_China_1969_7.25.2017.csv")) %>%
    dplyr::rename(
      gender = gender_it,
      age = age_it,
      height_dxa = height,
      weight_dxa = weight
    ) %>%
    mutate(
      aoi = android_fm / gynoid_fm,
      id = as.character(id),
      operator = NA,
      batch = NA
    ) %>%
    left_join(
      bc %>%
        dplyr::select(id, number, batch, operator)
    )

  taxa <- readRDS(file.path(data_dir, "taxa.rds")) %>%
    data.frame() %>%
    rownames_to_column("seq")
  taxa$seq_num <- colnames(seqtab)
  list(
    "taxa" = taxa,
    "seqtab" = seqtab,
    "bc" = bc,
    "bc_full" = bc_full,
    "tree" = readRDS(file.path(data_dir, "phylo_tree.rds"))
  )
}

#' Simulated data
#'
#' This creates the made-up data set, used just for others to experiment with
#' our code, without having to release the original (private) data.
#'
#' @param raw The output of the read_data() when given real data.
#' @return sim A dataset with the same structure as `raw`, but with all
#'   non-taxonomic data simulated.
create_sim <- function(raw) {
  sim <- list()

  # simulate the rsv matrix
  J <- 200
  sim$seqtab <- raw$seqtab[, 1:J]
  N <- nrow(sim$seqtab)
  rownames(sim$seqtab) <- paste0("sa", seq_len(N))
  for (j in seq_len(J)) {
    sim$seqtab[, j] <- rnbinom(N, 0.5, mu=15)
  }

  # simulate the body compositions
  sim$bc <- raw$bc
  sim$bc[, "id"] <- as.character(seq_len(nrow(sim$bc)))
  sim$bc[, "number"] <- as.character(seq_len(nrow(sim$bc)))
  sim$bc[, "gender"] <- sample(c("Male", "Female"), nrow(sim$bc), replace =  TRUE)
  sim$bc[, "age"] <- rpois(nrow(sim$bc), 50)
  sim$bc[, "batch"] <- 1
  sim$bc[, "operator"] <- 1
  for (j in 5:36) {
    sim$bc[, j] <- rnorm(nrow(sim$bc), 100, 20)
  }

  sim$bc_full <- sim$bc
  sim$tree <- raw$tree
  sim$taxa <- raw$taxa

  sim
}


#' Prepare sample data
#'
#' Reorders columns and converts to ratios.
#'
#' @param survey [data.frame] The raw body composition data frame.
#' @param scale_sample_data [bool] If TRUE, numeric columns will be centered and
#'   scaled (separately for each gender).
#' @return The body composition data with log transformed diet and ratio
#'   transformed lean and fat mass variables.
prep_sample <- function(survey, scale_sample_data = FALSE) {
  sample <- data.frame(survey) %>%
    mutate(
      id = as.character(id),
      gender = as.factor(gender)
    ) %>%
    select( ## reorder columns
      id, age, gender, height_dxa, weight_dxa,
      bmi, aoi, ends_with("_fm"), ends_with("_lm"),
      batch, operator
    )%>%
    mutate_at( ## fm ratios
      vars(ends_with("fm"), -total_fm),
      .funs = funs(. / total_fm)
    ) %>%
    mutate_at( ## lm ratios
      vars(ends_with("lm"), -total_lm),
      .funs = funs(. / total_lm)
    ) %>%
    mutate_at(
      vars(starts_with("diet_")),
      .funs = log
    ) %>%
    mutate( ## new variables
      fat_lean_ratio = total_fm / total_lm,
      total_fm = total_fm / (weight_dxa * 1000),
      total_lm = total_lm / (weight_dxa * 1000),
      batch = as.factor(batch),
      operator = as.factor(operator)
    )

  ## scale body measurements by gender
  if (scale_sample_data) {
    sample <- sample %>%
      group_by(gender) %>%
      mutate_if(is.numeric, safe_scale) %>%
      ungroup()
  }

  as.data.frame(sample)
}

#' Preprocess Taxa Table
#'
#' @param taxa [data.frame] The raw taxonomic data as a data.frame
#' @return Transformed version of the taxonomic table, includign a new "family"
#'   column with lumped factors, along with a "seq_num" column giving a concise
#'   species name for each sequence.
prepare_taxa <- function(taxa) {
  taxa <- taxa %>%
    mutate(
      family = fct_lump(Family, n = 9, ties.method = "first")
    )
  taxa[is.na(taxa[, "family"]), "family"] <- "Other"
  taxa$family <- factor(
    taxa$family,
    levels = names(sort(table(taxa$family), decreasing = TRUE))
  )
  taxa <- tax_table(as.matrix(taxa))
  taxa_names(taxa) <- taxa[, "seq_num"]
  taxa
}

#' Scaling that ignores NAs
#'
#' Scale a matrix ignoring any NA entries.
#'
#' @param x [matrix] The matrix to scale
#' @return x [matrix] The matrix with centered and scaled columns, ignoring any
#'   NA entries.
safe_scale <- function(x) {
  x[!is.finite(x)] <- NA
  z <- x - mean(x, na.rm = TRUE)
  z / sd(x, na.rm = TRUE)
}

#' VST transform a phyloseq object
vst_ps <- function(ps, sf_quantile = 0.95, ...) {
  if (sample_data(ps)$gender[1] == "Female") {
    fmla <- formula(~ gynoid_fm)
  } else {
    fmla <- formula(~ trunk_fm)
  }

  dds <- DESeqDataSetFromMatrix(
    countData = t(get_taxa(ps)),
    colData = sample_data(ps),
    design = fmla
  )

  qs <- apply(counts(dds), 2, quantile, sf_quantile)
  sizeFactors(dds) <- qs / exp(mean(log(qs)))
  varianceStabilizingTransformation(dds, ...)
}

#' Preprocess Microbiome + Body Composition Data
#'
#' Prepare the output of read_data() for use in actual analysis. This basically
#' wraps prepare_sample(), prepare_taxa() and some DESeq functions.
#'
#' @param seqtab [matrix] A matrix containing raw species counts across species.
#' @param bc [data.frame] The body composition data on the subset of people in
#'   the microbiome study.
#' @param bc_full [data.frame] The full body composition data, across all
#'   participants.
#' @param taxa [data.frame] The original taxa table for the species in seqtab.
#' @param tree [phylo] A phylogenetic tree relating species.
#' @param opts [list] A list giving preprocessing options. Supported options
#'   are,
#'     - scale_sample_data: Should body composition data be centered and scaled?
#'     - gender: If "Male" / "Female", only return male / female participants.
#'        Otherwise returns all.
#'     - filter_k: The proportion [in 0 to 1] of samples with more than a
#'         samples, needed for filtering.
#'     - filter_a: The threshold value used by the k over a filter.
#'     - outdir: To what directory should any intermediate rlog transformed data
#'         be saved?
#' @return resutls [list] A list containing preprocessed data. Includes the
#'   following elements,
#'     - ps [phyloseq]: A phyloseq object containing the merged otu, taxa, tree,
#'        and sample data.
#'     - bc_full [data.frame]: A data.frame containing body composition data for
#'        all study participants, not just those with microbiome data.
#'     - mseqtab [data.frame]: A "melted" version of the phyloseq object. Mainly
#'        useful for annotated plots.
process_data <- function(seqtab, bc, bc_full, taxa, tree, opts = list()) {
  ## preparing taxa and survey data
  opts <- merge_process_opts(opts)
  taxa <- prepare_taxa(taxa)
  taxa_names(tree) <- taxa[, "seq_num"]

  bc_full <- prep_sample(bc_full, opts$scale_sample_data) %>%
    left_join(bc %>% select(id, number))
  bc <- prep_sample(bc, opts$scale_sample_data) %>%
    left_join(bc %>% select(id, number))
  rownames(bc) <- bc$number

  ## combine into ps
  ps <- phyloseq(
    otu_table(seqtab, taxa_are_rows = FALSE),
    sample_data(bc),
    taxa,
    phy_tree(tree)
  )

  if (opts$gender == "Male") {
    ps <- ps %>%
      subset_samples(gender == "Male")
  } else if (opts$gender == "Female") {
    ps <- ps %>%
      subset_samples(gender == "Female")
  }

  ps <- ps %>%
    filter_taxa(
      function(x) {
        mean(x > opts$filter_a) > opts$filter_k
      },
      TRUE
    ) %>%
    transform_sample_counts(opts$transform_fun)

  if (opts$vst) {
    vst_dds <- do.call(vst_ps, c(ps, opts$vst_opts))
    otu_table(ps)@.Data <- t(assay(vst_dds))
  }

  if (opts$remove_na) {
    na_rows <- apply(sample_data(ps), 1, function(x) any(is.na(x)))
    ps <- ps %>%
      subset_samples(!na_rows)
  }

  list(
    "ps" = ps,
    "bc_full" = bc_full,
    "mseqtab" = melt_ps(ps)
  )
}

#' Melt a Phyloseq
#'
#' Transform a phyloseq object to a data.frame
#'
#' @param ps [phyloseq] A phyloseq object with counts, sample data, and
#'   taxonomic annotation to merge into a single data.frame
#' @return mps [data.frame] The melted phyloseq object as a single data.frame,
#'   with rows indeced by species x sample pairs.
melt_ps <- function(ps) {
  get_taxa(ps) %>%
    melt(varnames = c("number", "seq_num")) %>%
    left_join(sample_data(ps)) %>%
    left_join(data.frame(tax_table(ps)))
}

family_means <- function(mseqtab) {
  processed$mseqtab %>%
    group_by(family, number) %>%
    dplyr::summarise(family_mean = mean(value)) %>%
    spread(family, family_mean) %>%
    group_by(number) %>%
    dplyr::summarise(rl_ratio = tanh(Bacteroidaceae - Ruminococcaceae))
}
