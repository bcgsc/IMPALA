#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla
.libPaths("/projects/vporter_prj/R/x86_64-centos7-linux-gnu-library/4.0")

## ---------------------------------------------------------------------------
## Allelic Imbalance in Expression using MBASED, part 1
## Vanessa Porter, Oct. 2021
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(prob))
suppressMessages(library(tidyr))
suppressMessages(library(MBASED))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(stats))
suppressMessages(library(tibble))

## ---------------------------------------------------------------------------
## LOAD INPUT 
## ---------------------------------------------------------------------------

# Make help options
option_list = list(
  make_option(c("-p", "--phase"), type="character", default=NULL,
              help="Phased VCF file (from WhatsHap)", metavar="character"),
  make_option(c("-r", "--rna"), type="character", default=NULL,
              help="Tumour RNA vcf file (from Strelka2)", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default = NULL,
              help="Sample name from the RPKM matrix (HTMCP written like e.g. HTMCP.03.06.02109)", metavar="character"),
  make_option(c("-k", "--rpkm"), type="character", default = NULL,
              help="RPKM matrix", metavar="character"),
  make_option(c("-m", "--min"), type="numeric", default = 1,
              help="Minimum RPKM value", metavar="numeric"),
  make_option(c("-o", "--outdir"), type="character", default = "mbased_output",
              help="Output directory name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir
sample <- opt$sample
rpkm <- read.delim(opt$rpkm, header = T, stringsAsFactors = F) 
min <- opt$min

dir.create(out)

## ---------------------------------------------------------------------------
## USER FUNCTIONS
## ---------------------------------------------------------------------------

# extract info from a list
list_n_item <- function(list, n){
  sapply(list, `[`, n)
}

# define function to print out the summary of ASE results
summarizeASEResults_1s <- function(MBASEDOutput) {
  
  geneOutputDF <- data.frame(
    majorAlleleFrequency = assays(MBASEDOutput)$majorAlleleFrequency[,1],
    pValueASE = assays(MBASEDOutput)$pValueASE[,1],
    pValueHeterogeneity = assays(MBASEDOutput)$pValueHeterogeneity[,1])
  
  geneAllele <- as.data.frame(assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor) %>%
    rownames_to_column(var = "rowname") %>%
    dplyr::mutate(gene = unlist(lapply(strsplit(rowname, split = ":"),function(x){x = x[1]}))) %>%
    dplyr::group_by(gene) %>%
    summarise(allele1IsMajor = unique(mySample))
  
  geneOutputDF$allele1IsMajor <- geneAllele$allele1IsMajor[match(rownames(geneOutputDF), geneAllele$gene)]
  
  lociOutputGR <- rowRanges(metadata(MBASEDOutput)$locusSpecificResults)
  lociOutputGR$allele1IsMajor <- assays(metadata(MBASEDOutput)$locusSpecificResults)$allele1IsMajor[,1]
  lociOutputGR$MAF <- assays(metadata(MBASEDOutput)$locusSpecificResults)$MAF[,1]
  lociOutputList <- split(lociOutputGR, factor(lociOutputGR$aseID, levels=unique(lociOutputGR$aseID)))
  
  return(
    list(
      geneOutput=geneOutputDF,
      locusOutput=lociOutputList
    )
  )
}

## ---------------------------------------------------------------------------
## READ IN THE RNA SNV CALLS
## ---------------------------------------------------------------------------

# read in the RNA calls
rna <- read.delim("/projects/vporter_prj/tools/vporter-allelespecificexpression/output/HTMCP.03.06.02058/rna.isec.dna.snps.genes.vcf.gz", header = F, comment.char = "#", stringsAsFactors = F)
rna <- read.delim(opt$rna, header = F, comment.char = "#", stringsAsFactors = F)
colnames(rna) <- c("CHROM", "POS", "ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE1", "gene", "gene_biotype", "gene_locus") 
rna$variant <- paste0(rna$CHROM, ":", rna$POS)

# remove variants non-overlapping with genes of interest
rna_filt <- rna[rna$gene != ".",]

## ---------------------------------------------------------------------------
## EXTRACT REF/ALT READ COUNTS
## ---------------------------------------------------------------------------

# Extract and add the read counts
info <- strsplit(rna_filt$SAMPLE1, ":")
expr <- list_n_item(info, 6)
expr <- strsplit(expr, ",")
rna_filt$REF.COUNTS <- as.numeric(list_n_item(expr, 1))
rna_filt$ALT.COUNTS <- as.numeric(list_n_item(expr, 2))

## ---------------------------------------------------------------------------
## MBASED WITH OR WITHOUT PHASING
## ---------------------------------------------------------------------------

### WITH PHASING
if (!is.null(opt$phase)){
  
  ### 
  ### PHASING
  ###
  
  # WhatsHap phased VCF from ONT sequencing pipeline
  wh <- read.delim(opt$phase, header = F, comment.char = "#", stringsAsFactors = F)
  colnames(wh) <- c("CHROM",  "POS",  "ID",  "REF", "ALT", "QUAL", "FILTER",  "INFO", "FORMAT", "SAMPLE")
  wh$variant <- paste0(wh$CHROM, ":", wh$POS)
  
  # remove unphased variants - phased variants have the pipe "|" symbol in column 10 - and remove indels
  wh <- wh[grep("|", wh$SAMPLE, fixed=TRUE),]
  wh <- wh %>% dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1)
  
  # add genotype from the SAMPLE column as a new column
  info2 <- strsplit(wh$SAMPLE, ":")
  wh$GT <- list_n_item(info2, 1)
  
  # Add the genotype from WhatsHap 
  rna_filt$GT <- wh$GT[match(rna_filt$variant, wh$variant)]
  
  # annotate the phased variants as alleleA and alleleB
  rna_filt$alleleA <- ifelse(rna_filt$GT == "1|0", rna_filt$ALT, rna_filt$REF)
  rna_filt$alleleB <- ifelse(rna_filt$GT == "1|0", rna_filt$REF, rna_filt$ALT)
  
  # add the phased COUNTS variants as alleleA and alleleB
  rna_filt$alleleA.counts <- ifelse(rna_filt$GT == "1|0", rna_filt$ALT.COUNTS, rna_filt$REF.COUNTS)
  rna_filt$alleleB.counts <- ifelse(rna_filt$GT == "1|0", rna_filt$REF.COUNTS, rna_filt$ALT.COUNTS)
  
  # phased only variants
  rna_phased <- rna_filt[complete.cases(rna_filt),]
  
  # make SNV IDs
  rna_phased <- rna_phased %>%
    arrange(CHROM, POS) %>%
    group_by(gene) %>%
    mutate(label = paste0("SNV",1:n()))
  rna_phased$SNV.ID <- paste0(rna_phased$gene, ":", rna_phased$label)
  
  ### 
  ### MBASED
  ###
  
  print("Beginning MBASED ...")
  
  # make the GRanges object of the loci
  mySNVs <- GRanges(seqnames=rna_phased$CHROM,
                     ranges=IRanges(start=rna_phased$POS, width=1),
                     aseID=rna_phased$gene,
                     allele1=rna_phased$REF,
                     allele2=rna_phased$ALT)
  names(mySNVs) <- rna_phased$SNV.ID
  
  # create input RangedSummarizedExperiment object
  mySample <- SummarizedExperiment(
    assays=list(lociAllele1Counts=matrix(rna_phased$alleleA.counts,
                                         ncol=1,
                                         dimnames=list(names(mySNVs),'mySample')),
                lociAllele2Counts=matrix(rna_phased$alleleB.counts,
                                         ncol=1,
                                         dimnames=list(names(mySNVs),'mySample'))),
    rowRanges=mySNVs
  )
  
  # run MBASED
  ASEresults_1s_haplotypesKnown <- runMBASED(ASESummarizedExperiment=mySample,
                                             isPhased=TRUE,
                                             numSim=10^6,
                                             BPPARAM = SerialParam())
  
  saveRDS(ASEresults_1s_haplotypesKnown, file=paste0(out, "/ASEresults_1s_haplotypesKnown.rds"))
  
  # extract results
  results <- summarizeASEResults_1s(ASEresults_1s_haplotypesKnown)
  
  # adjust the pvalue with BH correction
  results$geneOutput$padj <- p.adjust(p = results$geneOutput$pValueASE, method = "BH")
  results$geneOutput$significance <- as.factor(ifelse(results$geneOutput$padj < 0.05, "padj < 0.05", "padj > 0.05"))
  results$geneOutput$gene <- rownames(results$geneOutput)

### WITHOUT PHASING
} else {
  
  # make SNV labels
  rna_filt <- rna_filt %>%
    arrange(CHROM, POS) %>%
    group_by(gene) %>%
    mutate(label = paste0("SNV",1:n()))
  rna_filt$SNV.ID <- paste0(rna_filt$gene, ":", rna_filt$label)
  
  ### 
  ### MBASED
  ###
  
  print("Beginning MBASED ...")
  
  # make the GRanges object of the loci
  mySNVs <- GRanges(seqnames=rna_filt$CHROM,
                    ranges=IRanges(start=rna_filt$POS, width=1),
                    aseID=rna_filt$gene,
                    allele1=rna_filt$REF,
                    allele2=rna_filt$ALT)
  names(mySNVs) <- rna_filt$SNV.ID
  
  ## create input RangedSummarizedExperiment object
  mySample <- SummarizedExperiment(
    assays=list(lociAllele1Counts=matrix(rna_filt$REF.COUNTS,
                                         ncol=1,
                                         dimnames=list(names(mySNVs),'mySample')),
                lociAllele2Counts=matrix(rna_filt$ALT.COUNTS,
                                         ncol=1,
                                         dimnames=list(names(mySNVs),'mySample'))),
    rowRanges=mySNVs
  )
  
  # run MBASED
  ASEresults_1s_haplotypesUnknown <- runMBASED(ASESummarizedExperiment=mySample,
                                               isPhased=FALSE,
                                               numSim=10^6,
                                               BPPARAM = SerialParam())
  saveRDS(ASEresults_1s_haplotypesUnknown, file=paste0(out, "/ASEresults_1s_haplotypesUnknown.rds"))
  
  # extract results
  results <- summarizeASEResults_1s(ASEresults_1s_haplotypesUnknown)
  
  # adjust the pvalue with BH correction
  results$geneOutput$padj <- p.adjust(p = results$geneOutput$pValueASE, method = "BH")
  results$geneOutput$significance <- as.factor(ifelse(results$geneOutput$padj < 0.05, "padj < 0.05", "padj > 0.05"))
  results$geneOutput$gene <- rownames(results$geneOutput)
  
} 

# save the results 
saveRDS(results, file=paste0(out, "/MBASEDresults.rds"))
print("Finished MBASED, adding expression ...")

## ---------------------------------------------------------------------------
## ADD THE RPKM AND FILTER THE ASE GENES
## ---------------------------------------------------------------------------

# add the RPKM of this sample
rpkm_sample <- rpkm[,c("gene", sample)] 

# expressed genes in the sample
results$geneOutput$RPKM <- rpkm_sample[match(results$geneOutput$gene, rpkm_sample$gene), 2]
results_filt <- results$geneOutput[results$geneOutput$RPKM > min, ]

# filter for genes that have an RPKM calculated
results_filt <- results_filt[complete.cases(results_filt),]

# MAF filter
results_filt$MAF <- as.factor(ifelse(results_filt$majorAlleleFrequency > 0.75, "MAF > 0.75", "MAF < 0.75"))

# add the locus
results_filt$gene_band <- rna_filt$gene_locus[match(results_filt$gene, rna_filt$gene)]
results_filt$gene_biotype <- rna_filt$gene_biotype[match(results_filt$gene, rna_filt$gene)]

# rearrange columns to a logical orger
results_filt <- results_filt[,c("gene", "gene_biotype", "gene_band", "RPKM", "allele1IsMajor","majorAlleleFrequency", 
                                "pValueASE", "pValueHeterogeneity", "padj",
                                "significance", "MAF")]

# save the data frame as a table
write.table(results_filt, paste0(out, "/MBASED_expr_gene_results.txt"), quote = F, col.names = T, row.names = F, sep = "\t")

print("Finished MBASED script")
