suppressMessages(library(optparse))
suppressMessages(library(tidyverse))

option_list = list(
  make_option(c("-e", "--expression_matrix"), type="character", default=NULL,
              help="Expression Matrix", metavar="character"),
  make_option(c("-t", "--tf"), type="character", default = NULL,
              help="Transcription Factor", metavar="character"),
  make_option(c("-a", "--allele1"), type="character", default = NULL,
              help="Allele 1", metavar="character"),
  make_option(c("-b", "--allele2"), type="character", default = NULL,
              help="Allele 2", metavar="character"),
  make_option(c("-i", "--id2gene"), type="character", default = NULL,
              help="id to gene tsv", metavar="character"),
  make_option(c("-m", "--min"), type="integer", default = 10,
              help="min expression level", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = NULL,
              help="output directory", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default = NULL,
              help="sample name", metavar="character")
)

# load in options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
min <- opt$min
outdir <- opt$outdir
sample <- opt$sample

## GET EXPRESSED TF
expressionMatrix <- read.delim(opt$expression_matrix, sep = "\t", header = T)
transcriptionFactor <- read.delim(opt$tf, sep = "\t", header = T) %>%
  separate(Transcription.factor, sep = ":", into=c("human", "transcriptionFactor")) %>%
  dplyr::select(-human)

expressionMatrix <- expressionMatrix[,c("gene", sample)]
colnames(expressionMatrix) <- c("gene", "expression")

expressionMatrix_filt <- expressionMatrix %>%
  dplyr::filter(gene %in% transcriptionFactor$transcriptionFactor) %>%
  dplyr::filter(expression > min)

tf_expression <- transcriptionFactor %>%
  left_join(expressionMatrix_filt, by=c("transcriptionFactor"="gene")) %>%
  dplyr::filter(!is.na(expression))

# GET ALLELE 1 TFBS
allele1MEME <- read.delim(opt$allele1) %>%
  dplyr::mutate(id = paste0(motif_id, "::", sequence_name)) %>%
  pull(id)

# GET ALLELE 2 TFBS
allele2MEME <- read.delim(opt$allele2) %>%
  dplyr::mutate(id = paste0(motif_id, "::", sequence_name)) %>%
  pull(id)

# GET GAIN OR LOSS OF TFBS
Gene <- read.delim(opt$id2gene, header = T)

difference <- rbind(tibble(id = setdiff(allele1MEME, allele2MEME), allele = 1),
                    tibble(id = setdiff(allele2MEME, allele1MEME), allele = 2)) %>%
  separate(id, sep = "::", into = c("motif_id", "sequence_id")) %>%
  left_join(Gene) %>%
  dplyr::filter(motif_id %in% tf_expression$Model) %>%
  left_join(tf_expression, by = c("motif_id"="Model")) %>%
  dplyr::select(transcriptionFactor, gene, allele) %>%
  write_tsv(paste0(outdir, "/motifDiff.tsv"))
