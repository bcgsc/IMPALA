suppressMessages(library(optparse))
suppressMessages(library(tidyverse))

option_list = list(
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help="Annotation file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = NULL,
              help="output directory", metavar="character")
)

# load in options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
annotation_path <- opt$annotation
out <- opt$out

read.delim(annotation_path, sep = "\t", header = T) %>%
  `colnames<-`(c("gene", "effect", "genotype", "filter" )) %>%
  dplyr::filter(filter == "PASS") %>%
  dplyr::filter(genotype == c("1|0", "0|1")) %>%
  dplyr::filter(effect != "stop_retained_variant") %>%
  dplyr::filter(grepl("stop", effect)) %>%
  distinct(gene, genotype, .keep_all = TRUE) %>%
  dplyr::mutate(stop_variant_allele = case_when(
    genotype == "0|1" ~ 2,
    TRUE ~ 1
  )) %>%
  dplyr::select(gene, effect, stop_variant_allele) %>%
  write.table(paste0(out, "/stop_variant.tsv"), 
              sep = "\t", quote = F, row.names = F, col.names = T)
  
  
  
