#!/gsc/software/linux-x86_64-centos7/R-4.1.3/bin/Rscript --vanilla
.libPaths("/home/glchang/R/x86_64-pc-linux-gnu-library/4.1")

suppressMessages(library("tidyverse"))
suppressMessages(library("optparse"))

option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="Sample name", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Gene expression data from RSEM", metavar="character"),
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help="Annotation file for ensembl id to hgnc name", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default =NULL,
              help="Output directory name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir
input <- opt$input
sample <- opt$sample
annotation <- read_delim(opt$annotation, col_names = F, show_col_types = F) %>%
  dplyr::select(X9, X4)

val <- read_delim(input, show_col_types = F) %>%
  dplyr::select(gene_id, TPM) %>%
  left_join(annotation, by = c("gene_id" = "X9")) %>%
  dplyr::select(X4, TPM) %>%
  `colnames<-`(c("gene", sample))

write.table(val, paste0(out, "/expression_matrix.tsv"), 
            col.names = T, row.names = F, quote = F, sep = "\t")
