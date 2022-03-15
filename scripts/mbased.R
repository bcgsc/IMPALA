#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla
.libPaths("/projects/vporter_prj/R/x86_64-centos7-linux-gnu-library/4.0")

## ---------------------------------------------------------------------------
## Allelic Imbalance in Expression using MBASED
## Vanessa Porter, Oct. 2021
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(prob))
suppressMessages(library(ggrepel))
suppressMessages(library(ggsci))
suppressMessages(library(tidyr))
suppressMessages(library(MBASED))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(stats))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tibble))
suppressMessages(library(chromPlot))
suppressMessages(library(networkD3))
suppressMessages(library(htmlwidgets))

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
print("Finished MBASED, beginning figures ...")

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
results_filt$gene_band <- all_genes$V8[match(results_filt$gene, all_genes$V4)]
results_filt$gene_biotype <- all_genes$V7[match(results_filt$gene, all_genes$V4)]

# rearrange columns to a logical orger
results_filt <- results_filt[,c("gene", "gene_biotype", "gene_band", "RPKM", "allele1IsMajor","majorAlleleFrequency", 
                                "pValueASE", "pValueHeterogeneity", "padj",
                                "significance", "MAF")]

# save the data frame as a table
write.table(results_filt, paste0(out, "/MBASED_expr_gene_results.txt"), quote = F, col.names = T, row.names = F, sep = "\t")

## ---------------------------------------------------------------------------
## FIGURES
## ---------------------------------------------------------------------------

### 
### SET UP THE DATAFRAME
###

df <- results_filt

# make a colour filter
df$colour_filt <- ifelse(df$padj < 0.05 & df$majorAlleleFrequency > 0.75, "MAF > 0.75 & padj < 0.05",
                         ifelse(df$padj > 0.05 & df$majorAlleleFrequency > 0.75, "MAF > 0.75 & padj > 0.05",
                                ifelse(df$padj > 0.05 & df$majorAlleleFrequency < 0.75, "MAF < 0.75 & padj > 0.05", "MAF < 0.75 & padj < 0.05")))

# add the chromosome
df$chr <- all_genes$V1[match(df$gene, all_genes$V4)]

# set the factor levels
df$colour_filt <- factor(df$colour_filt, levels = c("MAF > 0.75 & padj < 0.05", "MAF < 0.75 & padj < 0.05", 
                                                            "MAF > 0.75 & padj > 0.05", "MAF < 0.75 & padj > 0.05"))
df$chr <- factor(df$chr, levels = c(paste0("chr", 1:22), "chrX"))
df <- df[complete.cases(df),]

#### 
#### DOTPLOT
####

dotplot <- ggplot(df, aes(x = majorAlleleFrequency, y = padj, colour = colour_filt)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c("#e74645", "black", "black", "grey")) +
  theme_bw() + 
  geom_hline(yintercept = 0.05, linetype = 2) +
  geom_vline(xintercept = 0.75, linetype = 2) +
  geom_text(aes(label = paste0(table(colour_filt)["MAF > 0.75 & padj < 0.05"], " ASE genes"), x = 0.9, y = 0.75), size = 4.5, colour = "#e74645") +
  labs(x = "major allele frequency", y = "adjusted pvalue", colour = NULL) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.text = element_text(size = 10, colour = "black"))

ggsave(filename = paste0(out, "/aseGenesDot.pdf"), plot = dotplot, width = 5, height = 4, units = "in")

####
#### BARPLOT
####

barplot <- ggplot(df, aes(x = chr, fill = colour_filt)) +
  geom_bar() +
  scale_fill_manual(values = rev(c("#e0f0ea","#574f7d", "#95adbe", "#e74645"))) +
  theme_bw() +
  labs(x = "chromosome", y = "number of genes", fill = "ASE results")+
  theme(axis.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

ggsave(filename = paste0(out, "/aseGenesBar.pdf"), plot = barplot, width = 12, height = 5, units = "in")

####
#### CHROMPLOT
####

genes1 <- df[df$colour_filt == "MAF > 0.75 & padj < 0.05","gene"]
genes2 <- df[df$colour_filt != "MAF > 0.75 & padj < 0.05","gene"]

bed1 <- data.frame(Chrom = all_genes$V1[all_genes$V4 %in% genes1], 
                   Start = all_genes$V2[all_genes$V4 %in% genes1],
                   End = all_genes$V3[all_genes$V4 %in% genes1])
bed2 <- data.frame(Chrom = all_genes$V1[all_genes$V4 %in% genes2], 
                   Start = all_genes$V2[all_genes$V4 %in% genes2],
                   End = all_genes$V3[all_genes$V4 %in% genes2])

pdf(paste0(out, "/chromPlot.pdf"), width = 9, height = 8)
chromPlot(gaps=hg_gap, annot1=bed1, annot2 = bed2, colAnnot1 = "#e74645", colAnnot2 = "grey",  bands=hg_cytoBandIdeo, bin = 1000000, chrSide=c(-1,1,1,1,1,1,1,1))
dev.off()

####
#### SANKEY PLOT
####

# set filters on the RPKM matrix
rpkm_sample$gene_biotype <- all_genes$V7[match(rpkm_sample$gene, all_genes$V4)]
rpkm_sample_filt1 <- rpkm_sample[rpkm_sample$gene_biotype %in% c("lincRNA", "miRNA", "protein_coding"),]
rpkm_sample_filt2 <- rpkm_sample_filt1[rpkm_sample_filt1$HTMCP.03.06.02109 > 1,] 

# get the input values for the plot
a <- nrow(rpkm_sample_filt1)
b <- c(nrow(rpkm_sample_filt2),nrow(rpkm_sample_filt1[rpkm_sample_filt1$HTMCP.03.06.02109 <= 1,] ))
c <- c(nrow(results_filt),nrow(rpkm_sample_filt2[!rpkm_sample_filt2$gene %in% results_filt$gene,]))
d <- c(nrow(results_filt[results_filt$padj < 0.05 & results_filt$majorAlleleFrequency > 0.75,]),
       sum(nrow(results_filt[results_filt$padj >= 0.05 & results_filt$majorAlleleFrequency <= 0.75,]),
       nrow(results_filt[results_filt$padj >= 0.05 & results_filt$majorAlleleFrequency > 0.75,]),
       nrow(results_filt[results_filt$padj < 0.05 & results_filt$majorAlleleFrequency <= 0.75,])))

# create a connection data frame
links <- data.frame(
  source=c(rep(paste0("All Genes (n=", a, ")"), 2), rep(paste0("RPKM > 1 (n=", b[1], ")"), 2), rep(paste0("Phased Genes (n=", c[1], ")"), 2)),
  target=c(paste0("RPKM > 1 (n=", b[1], ")"), paste0("RPKM <= 1 (n=", b[2], ")"), 
           paste0("Phased Genes (n=", c[1], ")"), paste0("Unhased Genes (n=", c[2], ")"), 
           paste0("ASE Genes (n=", d[1], ")"), paste0("Biallelic Genes (n=", d[2], ")")), 
  value=c(b, c, d)
)


# create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# Reformat the links
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
sankey <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, fontSize = 18)

saveWidget(sankey, file=paste0(out, "/sankeyPlot.html"))

print("Figures completed")
