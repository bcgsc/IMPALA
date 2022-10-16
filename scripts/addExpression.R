
## ---------------------------------------------------------------------------
## ADD THE RPKM AND FILTER THE ASE GENES
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))

## ---------------------------------------------------------------------------
## OPTIONS
## ---------------------------------------------------------------------------

# Make help options
option_list = list(
  make_option(c("-b", "--mbased"), type="character", default=NULL,
              help="mbased rds file", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default = NULL,
              help="Sample name from the RPKM matrix (HTMCP written like e.g. HTMCP.03.06.02109)", metavar="character"),
  make_option(c("-r", "--rpkm"), type="character", default = NULL,
              help="RPKM matrix", metavar="character"),
  make_option(c("-m", "--min"), type="numeric", default = 1,
              help="Minimum RPKM value", metavar="numeric"),
  make_option(c("-t", "--maf_threshold"), type="numeric", default = 0.60,
              help="Threshold for MAF to consider as ASE", metavar="numeric"),
  make_option(c("-o", "--outdir"), type="character", default = "mBASED",
              help="Output directory name", metavar="character")
)

## ---------------------------------------------------------------------------
## VARIABLES
## ---------------------------------------------------------------------------

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

out <- opt$outdir
sample <- opt$sample
rpkm <- read.delim(opt$rpkm, header = T, stringsAsFactors = F) 
results <- readRDS(opt$mbased)
min <- opt$min
maf_threshold <- opt$maf_threshold 

## ---------------------------------------------------------------------------
## VARIABLES
## ---------------------------------------------------------------------------

print("Adding expression")

# fix sample name
sample <- ifelse(length(grep("-", sample)) == 0, sample, gsub("-", ".", sample))

# select the RPKM of this sample
rpkm_sample <- rpkm[,c("gene", sample)] 

# expressed genes in the sample
results$geneOutput$RPKM <- rpkm_sample[match(results$geneOutput$gene, gsub(" ", "", rpkm_sample$gene, fixed = TRUE)), 2]
results_filt <- results$geneOutput[results$geneOutput$RPKM > min, ]

# filter for genes that have an RPKM calculated
results_filt <- results_filt[!is.na(results_filt$RPKM),] 

# MAF filter
results_filt$MAF <- as.factor(ifelse(results_filt$majorAlleleFrequency > maf_threshold, paste0("MAF > ", maf_threshold), paste0("MAF < ", maf_threshold),))
results_filt$aseResults <- as.factor(ifelse(results_filt$majorAlleleFrequency > maf_threshold & results_filt$padj < 0.05, "ASE", "BAE"))

# rearrange columns to a logical order
results_filt <- results_filt[,c("gene", "geneBiotype", "RPKM", "allele1IsMajor","majorAlleleFrequency", 
                                "pValueASE", "pValueHeterogeneity", "padj",
                                "significance", "MAF", "aseResults")]

# save the data frame as a table
write.table(results_filt, paste0(out, "/MBASED_expr_gene_results.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
