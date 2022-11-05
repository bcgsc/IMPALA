suppressMessages(library(tidyverse))
suppressMessages(library(optparse))


# Make help options
option_list = list(
  make_option(c("-c", "--cnv"), type="character", default=NULL,
              help="CNV bed file", metavar="character"),
  make_option(c("-m", "--methyl"), type="character", default=NULL,
              help="allelic methylation bed file", metavar="character"),
  make_option(c("-a", "--ase"), type="character", default=NULL,
              help="mbased ASE result", metavar="character"),
  make_option(c("-g", "--cancer"), type="character", default=NULL,
              help="cancer gene", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="Sample name", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = NULL,
              help="Output directory name", metavar="character"),
  make_option(c("-t", "--tfbs"), type="character", default = NULL,
              help="Transcription Factor binding site", metavar="character"),
  make_option(c("-p", "--stop"), type="character", default = NULL,
              help="Stop Variants", metavar="character"),
  make_option(c("-n", "--snv"), type="character", default = NULL,
              help="Somatic SNV", metavar="character"),
  make_option(c("-i", "--indel"), type="character", default = NULL,
              help="Somatic Indel", metavar="character")
)


# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir
sample <- opt$sample

cnv_path <- opt$cnv
methyl_path <- opt$methyl
tfbs_path <- opt$tfbs
stop_path <- opt$stop
snv_path <- opt$snv
indel_path <- opt$indel
  
cancer <- read.delim(opt$cancer, sep = "\t", header = F, comment.char = "#") %>%
  pull()

ase <- read.delim(opt$ase, header = T, comment.char = "#", stringsAsFactors = F)

##########
# CNV
##########

if (is.null(cnv_path) | cnv_path == "") {
  cnv <- data.frame(gene = ase$gene)
} else {
  cnv <- read.delim(cnv_path, header = F, comment.char = "#") %>%
    dplyr::select(V4, V13, V14, V16) %>%
    dplyr::mutate(V13 = as.numeric(V13)) %>%
    dplyr::mutate(V14 = as.numeric(V14)) %>%
    dplyr::mutate(V16 = as.numeric(V16)) %>%
    `colnames<-`(c("gene", "cnv.A", "cnv.B", "expectedMAF")) %>%
    group_by(gene) %>% 
    summarize(cnv.A = mean(cnv.A), cnv.B = mean(cnv.B), expectedMAF = mean(expectedMAF)) %>%
    dplyr::mutate(cnv_state = case_when(
      cnv.A == 0 | cnv.B == 0 ~ "LOH",
      abs(cnv.A - cnv.B) <= 1 ~ "balance",
      TRUE ~ "imbalance"
    )) 
}

##########
# DMR
##########

if ( is.null(methyl_path) | methyl_path == "") {
  dmr <- data.frame(gene = ase$gene)
} else {
  dmr <- read.delim(methyl_path, header = F, comment.char = "#", stringsAsFactors = F) %>%
    dplyr::filter(V10 != ".") %>%
    dplyr::select(V4, V13, V14, V15) %>%
    `colnames<-`(c("gene", "methyl.A", "methyl.B", "diff.Methyl")) %>%
    dplyr::mutate(diff.Methyl = as.numeric(diff.Methyl)) %>%
    group_by(gene) %>%
    summarize(minimum = min(diff.Methyl), maximum = max(diff.Methyl)) %>%
    dplyr::mutate(methyl_state = case_when(
      (minimum < 0) & (maximum < 0) ~ minimum,
      (minimum > 0) & (maximum > 0) ~ maximum,
      TRUE ~ ifelse((abs(minimum) > maximum), minimum, maximum)
    )) %>%
    dplyr::select(gene, methyl_state)
}

########################################
# Transcription Factor Binding Site
########################################

if (is.null(tfbs_path) | tfbs_path == "") {
  tfbs <- data.frame(gene = ase$gene)
} else {
  tfbs <- read.delim(tfbs_path, header = T, comment.char = "#") %>%
    `colnames<-`(c("transcriptionFactor", "gene", "tf_allele")) %>%
    dplyr::select(gene, transcriptionFactor, tf_allele) %>%
    group_by(gene) %>% 
    mutate(transcriptionFactor = paste(transcriptionFactor, collapse=",")) %>%
    mutate(tf_allele = paste(tf_allele, collapse=",")) %>%
    distinct()
}

##########################
# Stop Variants
##########################

if (is.null(stop_path) | stop_path == "") {
  stopVar <- data.frame(gene = ase$gene)
} else {
  stopVar <- read.delim(stop_path, header = T, comment.char = "#") %>%
    dplyr::select(gene, stop_variant_allele)
}

##########################
# Somatic SNV
##########################

if (is.null(snv_path) | snv_path == "") {
  snv <- data.frame(gene = ase$gene)
} else {
  snv_gene <- read.delim(snv_path, header = F, comment.char = "#") %>%
    pull(V4) %>%
    unique()
  snv <- data.frame(gene = ase$gene) %>%
    dplyr::mutate(somaticSNV = gene %in% snv_gene)
}

##########################
# Somatic Indel
##########################

if (is.null(indel_path) | indel_path == "") {
  indel <- data.frame(gene = ase$gene)
} else {
  indel_gene <- read.delim(indel_path, header = F, comment.char = "#") %>%
    pull(V4) %>%
    unique()
  indel <- data.frame(gene = ase$gene) %>%
    dplyr::mutate(somaticIndel = gene %in% indel_gene)
}

##########################
# Summary Table 
##########################

summary_table <- ase %>%
  dplyr::select(gene, RPKM, allele1IsMajor, majorAlleleFrequency, padj, aseResults) %>%
  left_join(cnv, by = "gene") %>%
  left_join(dmr, by = "gene") %>%
  left_join(tfbs, by = "gene") %>%
  left_join(stopVar, by = "gene") %>%
  left_join(snv, by = "gene") %>%
  left_join(indel, by = "gene") %>%
  dplyr::mutate(cancer_gene = gene %in% cancer) %>%
  dplyr::mutate(sample = sample) %>%
  write.table(paste0(out, "/summaryTable.tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = T)

