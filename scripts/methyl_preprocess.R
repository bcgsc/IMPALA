suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

# Make help options
option_list = list(
  make_option(c("-m", "--methyl"), type="character", default=NULL,
              help="Allelic methylation (from Nanomethphase)", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = "mBASED",
              help="Output directory name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir


dmr <- read.delim(opt$methyl, header = T, comment.char = "#", stringsAsFactors = F)

dmr %>%
  dplyr::mutate(chr = case_when(
    grepl("^chr", chr) ~ chr,
    TRUE ~ paste0("chr", chr)
  )) %>%
  dplyr::select(chr, start, end, meanMethy1, meanMethy2, diff.Methy) %>%
  dplyr::mutate(start = gsub(" ", "", format(start, scientific=F), fixed = TRUE)) %>%
  dplyr::mutate(end = gsub(" ", "", format(end, scientific=F), fixed = TRUE)) %>%
  dplyr::mutate(meanMethy1 = format(meanMethy1, scientific = FALSE)) %>%
  dplyr::mutate(meanMethy2 = format(meanMethy2, scientific = FALSE)) %>%
  dplyr::mutate(diff.Methy = format(diff.Methy, scientific = FALSE)) %>%
  write.table(paste0(out, "/methyl.bed"), 
              quote = F, sep = "\t", row.names = F, col.names = F)
  
