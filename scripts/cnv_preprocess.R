suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

# Make help options
option_list = list(
  make_option(c("-c", "--cnv"), type="character", default=NULL,
              help="CNV file (from Ploidetect)", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = "mBASED",
              help="Output directory name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir


cnv <- read.delim(opt$cnv, header = T, comment.char = "#", stringsAsFactors = F)

cnv %>%
  dplyr::mutate(type = case_when(
    zygosity == "HOM" ~ "LOH",
    A != B ~ "imbalance",
    A == B ~ "balance"
  )) %>%
  dplyr::select(chr, pos, end, A, B, type) %>%
  dplyr::mutate(chr = case_when(
    grepl("^chr", chr) ~ chr,
    TRUE ~ paste0("chr", chr)
  )) %>%
  dplyr::mutate(expectedMAF = pmax(A, B)/(A + B)) %>%
  dplyr::mutate(pos = gsub(" ", "", format(pos, scientific=F), fixed = TRUE)) %>%
  dplyr::mutate(end = gsub(" ", "", format(end, scientific=F), fixed = TRUE)) %>%
  write.table(paste0(out, "/cnv.bed"), 
              quote = F, sep = "\t", row.names = F, col.names = F)

