library(ggplot2)
library(stringr)
library(optparse)

option_list <- list(
  make_option("--indir", type="character"),
  make_option("--outdir", type="character"),
  make_option("--assembly", type="character", default="Primary")
)

opt <- parse_args(OptionParser(option_list=option_list))

args <- commandArgs(trailingOnly = TRUE)

predector_dir <- args[1]
outdir <- args[2]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)




# Effector processing function


make_effector_table <- function(signalp_file, effectorp_file, tmhmm_file, assembly_name) {


  cat("Reading:", assembly_name, "\n")


  signalp <- read.delim(signalp_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  effp    <- read.delim(effectorp_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  tmhmm   <- read.delim(tmhmm_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)


  signalp$has_SP <- signalp$prediction == "SP"
  tmhmm$has_TMHMM <- tmhmm$pred_hel > 0


  effp$effector_class <- effp$prediction
  effp$effector_class[effp$prediction == "Apoplastic effector"] <- "Apoplastic"
  effp$effector_class[effp$prediction == "Cytoplasmic effector"] <- "Cytoplasmic"
  effp$effector_class[effp$prediction == "Cytoplasmic/apoplastic effector"] <- "Dual"
  effp$effector_class[effp$prediction == "Apoplastic/cytoplasmic effector"] <- "Dual"
  effp$effector_class[effp$prediction == "Non-effector"] <- "Non-effector"


  x <- merge(signalp, effp, by = "name", all.x = TRUE, suffixes = c("_signalp", "_effectorp"))
  x <- merge(x, tmhmm, by = "name", all.x = TRUE)


  x$assembly <- assembly_name
  x$has_TMHMM[is.na(x$has_TMHMM)] <- FALSE
  x$pred_hel[is.na(x$pred_hel)] <- 0


  x$final_secreted_class <- "No signal peptide"
  x$final_secreted_class[x$has_SP & x$has_TMHMM] <- "SP + TMHMM excluded"
  x$final_secreted_class[x$has_SP & !x$has_TMHMM] <- x$effector_class[x$has_SP & !x$has_TMHMM]


  x$retained_for_effector_plot <- x$has_SP & !x$has_TMHMM


  return(x)
}


# Predector output files are our input files




effector_table <- make_effector_table(
  file.path(outdir, "Pri/S1W6_primary_proteins_clean-signalp6.tsv"),
  file.path(outdir, "Pri/S1W6_primary_proteins_clean-effectorp3.tsv"),
  file.path(outdir, "Pri/S1W6_primary_proteins_clean-tmhmm.tsv"),
  "Primary"
)


write.table(
  effector_table,
  file.path(outdir, "Effector_secreted_effector_ID_table.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)




# 1. Effector class summary for retained secreted and no-TM proteins


plot_table <- effector_table[effector_table$retained_for_effector_plot, ]


effector_counts <- as.data.frame(table(
  plot_table$assembly,
  plot_table$effector_class
))


colnames(effector_counts) <- c("assembly", "effector_class", "n")


effector_counts <- effector_counts[effector_counts$n > 0, ]


write.table(
  effector_counts,
  file.path(outdir, "Effector_stacked_bar_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


# 2. Effector class barplot




p_effector <- ggplot(
  effector_counts,
  aes(x = assembly, y = n, fill = effector_class)
) +
  geom_col(color = "black") +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Number of proteins",
    fill = "EffectorP class"
  )


ggsave(
  file.path(outdir, "Effector_stacked_bar_plot.pdf"),
  p_effector,
  width = 5,
  height = 4
)


ggsave(
  file.path(outdir, "Effector_stacked_bar_plot.png"),
  p_effector,
  width = 5,
  height = 4,
  dpi = 300
)




# 3. CAZyme class summary table




cazy_file <- file.path(outdir, "Pri/S1W6_primary_proteins_clean-dbcan.tsv")


if (!file.exists(cazy_file)) {
  cazy_candidates <- list.files(
    file.path(outdir, "Pri"),
    pattern = "dbcan|dbCAN|cazy|CAZy|hmmer|HMMER",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )


  if (length(cazy_candidates) == 0) {
    stop("No CAZyme/dbCAN/HMMER file found in Pri/.")
  }


  cazy_file <- cazy_candidates[1]
}


cat("Reading CAZyme file:\n", cazy_file, "\n")


cazy <- read.delim(cazy_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)


if (!"hmm" %in% colnames(cazy)) {
  stop("Expected a column named 'hmm'. Columns are: ", paste(colnames(cazy), collapse = ", "))
}


cazy$cazy_family <- as.character(cazy$hmm)
cazy$cazy_class <- str_extract(cazy$cazy_family, "^(AA|CBM|CE|GH|GT|PL)")


cazy_counts <- as.data.frame(table(cazy$cazy_class))
colnames(cazy_counts) <- c("cazy_class", "n")
cazy_counts <- cazy_counts[!is.na(cazy_counts$cazy_class) & cazy_counts$cazy_class != "", ]
cazy_counts$assembly <- "Primary"


cazy_counts <- cazy_counts[, c("assembly", "cazy_class", "n")]


write.table(
  cazy,
  file.path(outdir, "CAZyme_full_table.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


write.table(
  cazy_counts,
  file.path(outdir, "CAZyme_stacked_bar_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)




# 4. CAZyme class barplot


p_cazy <- ggplot(
  cazy_counts,
  aes(x = assembly, y = n, fill = cazy_class)
) +
  geom_col(color = "black") +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Number of proteins",
    fill = "CAZyme class"
  )


ggsave(
  file.path(outdir, "CAZyme_stacked_bar_plot.pdf"),
  p_cazy,
  width = 5,
  height = 4
)


ggsave(
  file.path(outdir, "CAZyme_stacked_bar_plot.png"),
  p_cazy,
  width = 5,
  height = 4,
  dpi = 300
)




# 5. Overall summary table




overall_summary <- data.frame(
  assembly = "Primary",
  total_genes = nrow(effector_table),
  number_SP = sum(effector_table$has_SP, na.rm = TRUE),
  non_secreted_number = sum(!effector_table$has_SP, na.rm = TRUE),
  number_with_TMHMM_domains = sum(effector_table$has_TMHMM, na.rm = TRUE),
  secreted_no_TMHMM = sum(effector_table$has_SP & !effector_table$has_TMHMM, na.rm = TRUE),
  candidate_effectors_apoplastic_plus_dual = sum(
    effector_table$retained_for_effector_plot &
      effector_table$effector_class %in% c("Apoplastic", "Dual"),
    na.rm = TRUE
  ),
  total_CAZyme_hits = nrow(cazy),
  total_unique_CAZyme_genes = length(unique(cazy$query))
)


write.table(
  overall_summary,
  file.path(outdir, "Effector_CAZyme_overall_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


# Print results


print(effector_counts)
print(cazy_counts)
print(overall_summary)


cat("\nDONE. Files written to:\n", outdir, "\n")