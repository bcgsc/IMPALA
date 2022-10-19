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
              help="Output directory name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir
sample <- opt$sample
cancer <- read.delim(opt$cancer, sep = "\t", header = F, comment.char = "#") %>%
  pull()

cnv <- read.delim(opt$cnv, header = F, comment.char = "#") %>%
  dplyr::select(V4, V12, V13, V15) %>%
  `colnames<-`(c("gene", "cnv.A", "cnv.B", "expectedMAF")) %>%
  group_by(gene) %>% 
  summarize(cnv.A = mean(cnv.A), cnv.B = mean(cnv.B), expectedMAF = mean(expectedMAF)) %>%
  dplyr::mutate(cnv_state = case_when(
    cnv.A == 0 | cnv.B == 0 ~ "LOH",
    abs(cnv.A - cnv.B) <= 1 ~ "balance",
    TRUE ~ "imbalance"
  )) 

dmr <- read.delim(opt$methyl, header = F, comment.char = "#", stringsAsFactors = F) %>%
  dplyr::select(V4, V12, V13, V14) %>%
  `colnames<-`(c("gene", "methyl.A", "methyl.B", "diff.Methyl")) %>%
  group_by(gene) %>%
  summarize(minimum = min(diff.Methyl), maximum = max(diff.Methyl)) %>%
  dplyr::mutate(methyl_state = case_when(
    (minimum < 0) & (maximum < 0) ~ minimum,
    (minimum > 0) & (maximum > 0) ~ maximum,
    TRUE ~ ifelse((abs(minimum) > maximum), minimum, maximum)
  )) %>%
  dplyr::select(gene, methyl_state)



ase <- read.delim(opt$ase, header = T, comment.char = "#", stringsAsFactors = F)

summary_table <- ase %>%
  dplyr::select(gene, RPKM, allele1IsMajor, majorAlleleFrequency, padj, aseResults) %>%
  left_join(cnv, by = "gene") %>%
  left_join(dmr, by = "gene") %>%
  dplyr::mutate(sample = sample) %>%
  dplyr::mutate(cancer_gene = gene %in% cancer) %>%
  write.table(paste0(out, "/summaryTable.tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = T)

