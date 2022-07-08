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
  make_option(c("-o", "--outdir"), type="character", default = "mBASED",
              help="Output directory name", metavar="character")
)

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir

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
rna_filt <- read.delim(opt$rna, header = T, comment.char = "#", stringsAsFactors = F)
colnames(rna_filt) <- c("CHROM", "POS", "AD","REF","ALT","gene", "gene_biotype") 
rna_filt$variant <- paste0(rna_filt$CHROM, ":", rna_filt$POS)


## ---------------------------------------------------------------------------
## EXTRACT REF/ALT READ COUNTS
## ---------------------------------------------------------------------------

# Extract and add the read counts
expr <- strsplit(rna_filt$AD, ",")
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
  
  # Find unphased genes with one variant (test)
  singleUnphased <- rna_filt %>%
    mutate(phase = variant %in% wh$variant) %>%
    left_join(rna_filt %>% group_by(gene) %>% summarize(n=n())) %>%
    dplyr::filter(!phase & n == 1)
  
  # Add genotype to unphased gene with one variant (test)
  rna_filt$GT[which(rna_filt$variant %in% singleUnphased$variant)] <- "1|0"
  
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
                                             BPPARAM = MulticoreParam(workers = 20))
  
  saveRDS(ASEresults_1s_haplotypesKnown, file=paste0(out, "/ASEresults_1s_haplotypesKnown.rds"))
  # extract results
  results <- summarizeASEResults_1s(ASEresults_1s_haplotypesKnown)
  
  # adjust the pvalue with BH correction
  results$geneOutput$padj <- p.adjust(p = results$geneOutput$pValueASE, method = "BH")
  results$geneOutput$significance <- as.factor(ifelse(results$geneOutput$padj < 0.05, "padj < 0.05", "padj > 0.05"))
  results$geneOutput$gene <- rownames(results$geneOutput)
  
  results$geneOutput$allele1IsMajor[results$geneOutput$gene %in% singleUnphased$gene] = NA
  
  # add the locus
  results$geneOutput$geneBiotype <- rna_filt$gene_biotype[match(results$geneOutput$gene, rna_filt$gene)]

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
                                               BPPARAM = MulticoreParam(workers = 20))
  saveRDS(ASEresults_1s_haplotypesUnknown, file=paste0(out, "/ASEresults_1s_haplotypesUnknown.rds"))
  
  # extract results
  results <- summarizeASEResults_1s(ASEresults_1s_haplotypesUnknown)
  
  # adjust the pvalue with BH correction
  results$geneOutput$padj <- p.adjust(p = results$geneOutput$pValueASE, method = "BH")
  results$geneOutput$significance <- as.factor(ifelse(results$geneOutput$padj < 0.05, "padj < 0.05", "padj > 0.05"))
  results$geneOutput$gene <- rownames(results$geneOutput)
  
  # add the locus
  results$geneOutput$geneBiotype <- rna_filt$gene_biotype[match(results$geneOutput$gene, rna_filt$gene)]
  
} 

# save the results 
saveRDS(results, file=paste0(out, "/MBASEDresults.rds"))
print("Finished MBASED")

