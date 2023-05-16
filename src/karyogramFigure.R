#!/gsc/software/linux-x86_64-centos7/R-4.0.2/bin/Rscript --vanilla
## ---------------------------------------------------------------------------
## GENERATING THE KARYOPLOT FROM ASE DATA
## Vanessa Porter, Dec. 2022
## ---------------------------------------------------------------------------
library(optparse)

## ---------------------------------------------------------------------------
## OPT PARSE
## ---------------------------------------------------------------------------

# options
option_list = list(
  make_option(c("-c", "--chromSize"), type="character",  
              help="chromosome sizes", metavar="character"),
  make_option(c("-p", "--centPos"), type="character", 
              help="centred centromere positions", metavar="character"),
  make_option(c("-n", "--cna"), type="character", 
              help="copy number segment file (preferably condensed)",  metavar="character"),
  make_option(c("-d", "--dmr"), type="character", 
              help="DMR positions", metavar="character"),
  make_option(c("-a", "--ase"), type="character",
              help="ASE summary table file", metavar="character"),
  make_option(c("-g", "--genes"), type="character",
              help="Gene annotation file", metavar="character"),
  make_option(c("-o", "--out"), type="character", 
              help="output file prefix", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Note: these packages need to be installed.
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(cowplot))

# chromosome size/position info
#chromSize <- read.delim("/projects/hpv_nanopore_prj/refs/hg38_no_alt_TCGA_HTMCP_HPVs_chromSizes.txt", header = F)
#centPos <-  read.delim("/projects/hpv_nanopore_prj/refs/hg38_centromere_positions_merged.bed", header = F)
#genes <- read.delim("/projects/hpv_nanopore_prj/htmcp/ase/pull_trial/vporter-allelespecificexpression/output/HTMCP.03.06.02058/3_cancer/raw/gene_annotation.bed", header = F)
chromSize <- read.delim(opt$chromSize, header = F)
centPos <-  read.delim(opt$centPos, header = F)
genes <- read.delim(opt$genes, header = F)

## ---------------------------------------------------------------------------
## POSITION CHROMOSOME INFO
## ---------------------------------------------------------------------------

# subset to the main chromosomes
chromSize <- chromSize[chromSize$V1 %in% c(paste0("chr", 1:22), "chrX"),]

# Rename columns to chromosome and size
colnames(chromSize) <- c("chr","size")

# Reorder levels for plotting
chromSize$chr <- factor(chromSize$chr,levels=c(paste0("chr", 1:22), "chrX"))
chromSize <- chromSize[chromSize$chr %in% c(paste0("chr", 1:22), "chrX"),]

# Divide by 1Mb to clean up axis
chromSize$size <- chromSize$size/1000000

# centromere mapping
colnames(centPos) <- c("chr", "start", "end")
centPos$chr <- factor(centPos$chr,levels=c(paste0("chr", 1:22), "chrX"))
centPos$centre <- centPos$start + ((centPos$end - centPos$start)/2)
centPos$centre <- centPos$centre/1000000

## ---------------------------------------------------------------------------
## COPY NUMBER
## ---------------------------------------------------------------------------
if (!is.null(opt$cna) & opt$cna != ""){
  #cna <- read.delim("/projects/hpv_nanopore_prj/htmcp/ploidetect/illumina/Ploidetect-pipeline/ploidetect_out/HTMCP-03-06-02058/A37261_A37189/cna_condensed.txt", header = T)
  cna <- read.delim(opt$cna, header = T)
  cna$chr <- paste0("chr", cna$chr)
  
  # rearrange to make a bed file
  cna_bed <- cna[,c("chr", "pos", "end", "CN", "zygosity", "A", "B")]
  
  # categorize copy number
  cna_bed <- cna_bed %>%
    mutate(CN.Status = 
             case_when(
               zygosity == "HOM" ~ "LOH",
               A > B ~ "imbalance",
               TRUE ~ "balance"
             ))
  
  # Divide by 1Mb for axis
  cna_bed$end <- cna_bed$end/1000000
  cna_bed$pos <- cna_bed$pos/1000000
  
  # Change to factor and reorder levels
  cna_bed$chr <- factor(cna_bed$chr,levels=c(paste0("chr", 1:22), "chrX"))
  
  cnaLOH <- cna_bed %>% filter(CN.Status == "LOH")
  cnaGAIN <- cna_bed %>% filter(CN.Status == "imbalance") 
}

## ---------------------------------------------------------------------------
## DIFFERENTIAL METHYLATION
## ---------------------------------------------------------------------------

if (!is.null(opt$dmr) & opt$dmr != ""){
  #dmr <- read.delim("/projects/hpv_nanopore_prj/htmcp/call_integration/output/HTMCP-03-06-02058/methylation/diff_meth.csv", header = T)
  dmr <- read.delim(opt$dmr, header = T)
  
  # Divide by 1Mb for axis
  dmr$start <- dmr$start/1000000
  dmr$end <- dmr$end/1000000
  dmr$middle <- (dmr$start + dmr$end) / 2
  
  # Change to factor and reorder levels
  dmr$chr <- factor(dmr$chr, levels=c(paste0("chr", 1:22), "chrX"))
  
  # count in 1Mb bins
  dmrCount <- data.frame(table(as.factor(paste0(dmr$chr, ":", as.integer(dmr$middle)))))
  
  # split the chromosome name and bin position
  dmrPlot <- separate(dmrCount, col = Var1, into = c("chr", "pos"), sep = ":", remove = T)
  
  # scale to fit the plot - i.e. make the maximum width 0.65
  maxDMR <- max(dmrPlot$Freq)
  dmrPlot$percMax <- dmrPlot$Freq/maxDMR
  dmrPlot$percMax <- dmrPlot$percMax * 0.65
  
  # Change to factor and reorder levels
  dmrPlot$chr <- factor(dmrPlot$chr, levels=c(paste0("chr", 1:22), "chrX"))
  dmrPlot$pos <- as.numeric(dmrPlot$pos)
}

## ---------------------------------------------------------------------------
## ASE GENE HISTOGRAM
## ---------------------------------------------------------------------------

#ase <- read.delim("/projects/hpv_nanopore_prj/htmcp/ase/pull_trial/vporter-allelespecificexpression/output/HTMCP.03.06.02058/summaryTable.tsv", header = T)
ase <- read.delim(opt$ase, header = T)

# filter for ASE genes
ase <- ase[ase$aseResults == "ASE",]

# get gene positions
ase$chr <- genes$V1[match(ase$gene, genes$V4)]
ase$start <- genes$V2[match(ase$gene, genes$V4)]
ase$end <- genes$V3[match(ase$gene, genes$V4)]

# get the middle of the gene for plotting
ase$middle <- (ase$start + ase$end)/2

# get the data frame ready for plotting 
ase <- ase[,c("chr", "start", "end","middle")]
ase <- ase[complete.cases(ase),]

# Divide by 1Mb for axis
ase$middle <- ase$middle/1000000

# count in 1Mb bins
aseCount <- data.frame(table(as.factor(paste0(ase$chr, ":", as.integer(ase$middle)))))

# split the chromosome name and bin position
asePlot <- separate(aseCount[aseCount$Var1 != "NA:NA",], col = Var1, into = c("chr", "pos"), sep = ":", remove = T)

# scale to fit the plot - i.e. make the maximum width 0.65
maxASE <- max(asePlot$Freq)
asePlot$percMax <- asePlot$Freq/maxASE
asePlot$percMax <- asePlot$percMax * 0.65

# Change to factor and reorder levels
asePlot$chr <- factor(asePlot$chr, levels=c(paste0("chr", 1:22), "chrX"))
asePlot$pos <- as.numeric(asePlot$pos)

## ---------------------------------------------------------------------------
## PLOT OPTIONS
## ---------------------------------------------------------------------------

##### CNV AND DMRs AVAILABLE

if (!is.null(opt$dmr) & opt$dmr != "" & !is.null(opt$cna) & opt$cna != ""){
  # legend
  adL <- data.frame(xmin = c(9.7, 9.7, 9.7, 10.3,9.7,10.3,7.1,7.1), xmax = c(10.35, 10.35,9.75,10.35,9.75,10.35,7.26,7.26), ymin = c(240,210,238,238,208,208,235,205), ymax = c(243,213,243,243,213,213,245,215),
                    fill = c("ase","dmr","ase","ase","dmr","dmr", "gain", "loh"))
  
  adW <- data.frame(x = c(11.5, 11.2,8.25,7.6,10,10), y = c(240,210,240,210,250,220),
                    label = c("ASE Gene Density","DMR Density", "Imbalanced CNV", "LOH", as.character(c(maxASE,maxDMR))))
  
  # plot
  # chromosomes 1 - 12
  p1 <- ggplot() +
    # chromosome bars
    geom_segment(data = chromSize %>% filter(chr %in% paste0("chr", 1:12)), aes(x = chr, xend = chr, y = 0, yend = size), 
                 lineend = "round", color = "lightgrey", size = 5) +
    # LOH
    geom_rect(data = cnaLOH %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = as.integer(chr) - 0.08, xmax = as.integer(chr) + 0.08, ymin = pos, ymax = end),
              fill="#94d2bd",size = 0.2) +
    # Imbalanced CNV
    geom_rect(data = cnaGAIN %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = as.integer(chr) - 0.08, xmax = as.integer(chr) + 0.08, ymin = pos, ymax = end),
              fill="#ee9b00",size = 0.2) +
    # ASE genes
    geom_rect(data = asePlot %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = as.integer(chr) + 0.1, xmax = (as.integer(chr) + 0.1 + percMax), ymin = pos, ymax = pos+1),
              fill = "#005f73", size = 0.25) +
    # DMRs
    geom_rect(data = dmrPlot %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = (as.integer(chr) - 0.1 - percMax), xmax = as.integer(chr) - 0.1, ymin = pos, ymax = pos+1),
              fill = "#ae2012", size = 0.25) +
    # centromeres
    geom_point(data = centPos %>% filter(chr %in% paste0("chr", 1:12)), aes(x = chr, y = centre), 
               size = 5, colour = "gray") +
    # legend bars
    geom_rect(data = adL, 
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              size = 0.25) +
    # legend text
    geom_text(data = adW, 
              aes(x = x, y = y, label = label))+
    scale_fill_manual(values = c("#005f73","#ae2012","#ee9b00","#94d2bd")) +
    ylim(0, 250) +
    theme_classic() +
    theme(text = element_text(size=15),axis.line=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")+
    labs(x=NULL,y="Chromosome Size (Mb)")
  
  # chromosomes 13 - 22 + X
  
  # very annoying but you have to filter all the dataframes or else the factor levels won't match the integer value
  chromSizeFilt <- chromSize %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  chromSizeFilt$chr <- factor(chromSizeFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  centPosFilt <- centPos %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  centPosFilt$chr <- factor(centPosFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  cnaLOHFilt <- cnaLOH %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  cnaLOHFilt$chr <- factor(cnaLOHFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  cnaGAINFilt <- cnaGAIN %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  cnaGAINFilt$chr <- factor(cnaGAINFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  asePlotFilt <- asePlot %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  asePlotFilt$chr <- factor(asePlotFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  dmrPlotFilt <- dmrPlot %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  dmrPlotFilt$chr <- factor(dmrPlotFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  
  # chromosomes 13-22+X
  p2 <- ggplot() +
    # chromosome bars
    geom_segment(data = chromSizeFilt, aes(x = chr, xend = chr, y = 0, yend = size), 
                 lineend = "round", color = "lightgrey", size = 5) +
    # LOH
    geom_rect(data = cnaLOHFilt, 
              aes(xmin = as.integer(chr) - 0.08, xmax = as.integer(chr) + 0.08, ymin = pos, ymax = end),
              fill="#94d2bd",size = 0.2) +
    # Imbalanced CNV
    geom_rect(data = cnaGAINFilt, 
              aes(xmin = as.integer(chr) - 0.08, xmax = as.integer(chr) + 0.08, ymin = pos, ymax = end),
              fill="#ee9b00",size = 0.2) +
    # ASE genes
    geom_rect(data = asePlotFilt, 
              aes(xmin = as.integer(chr) + 0.1, xmax = (as.integer(chr) + 0.1 + percMax), ymin = pos, ymax = pos+1),
              fill = "#005f73", size = 0.25) +
    # DMRs
    geom_rect(data = dmrPlotFilt, 
              aes(xmin = (as.integer(chr) - 0.1 - percMax), xmax = as.integer(chr) - 0.1, ymin = pos, ymax = pos+1),
              fill = "#ae2012", size = 0.25) +
    # centromeres
    geom_point(data = centPosFilt, aes(x = chr, y = centre), 
               size = 5, colour = "black") +
    ylim(0, 250) +
    theme_classic() +
    theme(text = element_text(size=15),axis.line=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")+
    labs(x=NULL,y=NULL)
  
} else if (!is.null(opt$cna) & opt$cna != ""){ ### CNVs BUT NO DMRs
  # legend
  adL <- data.frame(xmin = c(9.7, 9.7, 10.3,7.1,7.1), xmax = c(10.35, 9.75,10.35,7.26,7.26), ymin = c(240,238,238,235,205), ymax = c(243,243,243,245,215),
                    fill = c("ase","ase","ase","gain", "loh"))
  
  adW <- data.frame(x = c(11.5,8.25,7.6,10), y = c(240,240,210,250),
                    label = c("ASE Gene Density", "Imbalanced CNV", "LOH", as.character(c(maxASE))))
  
  # plot
  # chromosomes 1 - 12
  p1 <- ggplot() +
    # chromosome bars
    geom_segment(data = chromSize %>% filter(chr %in% paste0("chr", 1:12)), aes(x = chr, xend = chr, y = 0, yend = size), 
                 lineend = "round", color = "lightgrey", size = 5) +
    # LOH
    geom_rect(data = cnaLOH %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = as.integer(chr) - 0.08, xmax = as.integer(chr) + 0.08, ymin = pos, ymax = end),
              fill="#94d2bd",size = 0.2) +
    # Imbalanced CNV
    geom_rect(data = cnaGAIN %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = as.integer(chr) - 0.08, xmax = as.integer(chr) + 0.08, ymin = pos, ymax = end),
              fill="#ee9b00",size = 0.2) +
    # ASE genes
    geom_rect(data = asePlot %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = as.integer(chr) + 0.1, xmax = (as.integer(chr) + 0.1 + percMax), ymin = pos, ymax = pos+1),
              fill = "#005f73", size = 0.25) +
    # centromeres
    geom_point(data = centPos %>% filter(chr %in% paste0("chr", 1:12)), aes(x = chr, y = centre), 
               size = 5, colour = "gray") +
    # legend bars
    geom_rect(data = adL, 
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              size = 0.25) +
    # legend text
    geom_text(data = adW, 
              aes(x = x, y = y, label = label))+
    scale_fill_manual(values = c("#005f73","#ee9b00","#94d2bd")) +
    ylim(0, 250) +
    theme_classic() +
    theme(text = element_text(size=15),axis.line=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")+
    labs(x=NULL,y="Chromosome Size (Mb)")
  
  # chromosomes 13 - 22 + X
  
  # very annoying but you have to filter all the dataframes or else the factor levels won't match the integer value
  chromSizeFilt <- chromSize %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  chromSizeFilt$chr <- factor(chromSizeFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  centPosFilt <- centPos %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  centPosFilt$chr <- factor(centPosFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  cnaLOHFilt <- cnaLOH %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  cnaLOHFilt$chr <- factor(cnaLOHFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  cnaGAINFilt <- cnaGAIN %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  cnaGAINFilt$chr <- factor(cnaGAINFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  asePlotFilt <- asePlot %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  asePlotFilt$chr <- factor(asePlotFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  
  # chromosomes 13-22+X
  p2 <- ggplot() +
    # chromosome bars
    geom_segment(data = chromSizeFilt, aes(x = chr, xend = chr, y = 0, yend = size), 
                 lineend = "round", color = "lightgrey", size = 5) +
    # LOH
    geom_rect(data = cnaLOHFilt, 
              aes(xmin = as.integer(chr) - 0.08, xmax = as.integer(chr) + 0.08, ymin = pos, ymax = end),
              fill="#94d2bd",size = 0.2) +
    # Imbalanced CNV
    geom_rect(data = cnaGAINFilt, 
              aes(xmin = as.integer(chr) - 0.08, xmax = as.integer(chr) + 0.08, ymin = pos, ymax = end),
              fill="#ee9b00",size = 0.2) +
    # ASE genes
    geom_rect(data = asePlotFilt, 
              aes(xmin = as.integer(chr) + 0.1, xmax = (as.integer(chr) + 0.1 + percMax), ymin = pos, ymax = pos+1),
              fill = "#005f73", size = 0.25) +
    # centromeres
    geom_point(data = centPosFilt, aes(x = chr, y = centre), 
               size = 5, colour = "black") +
    ylim(0, 250) +
    theme_classic() +
    theme(text = element_text(size=15),axis.line=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")+
    labs(x=NULL,y=NULL)
  
} else if (!is.null(opt$dmr) & opt$dmr != ""){ ## DMRs BUT NO CNV
  # legend
  adL <- data.frame(xmin = c(9.7, 9.7, 9.7, 10.3,9.7,10.3), xmax = c(10.35, 10.35,9.75,10.35,9.75,10.35), ymin = c(240,210,238,238,208,208), ymax = c(243,213,243,243,213,213),
                    fill = c("ase","dmr","ase","ase","dmr","dmr"))
  
  adW <- data.frame(x = c(11.5, 11.2,10,10), y = c(240,210,250,220),
                    label = c("ASE Gene Density","DMR Density", as.character(c(maxASE,maxDMR))))
  
  # plot
  # chromosomes 1 - 12
  p1 <- ggplot() +
    # chromosome bars
    geom_segment(data = chromSize %>% filter(chr %in% paste0("chr", 1:12)), aes(x = chr, xend = chr, y = 0, yend = size), 
                 lineend = "round", color = "lightgrey", size = 5) +
    # ASE genes
    geom_rect(data = asePlot %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = as.integer(chr) + 0.1, xmax = (as.integer(chr) + 0.1 + percMax), ymin = pos, ymax = pos+1),
              fill = "#005f73", size = 0.25) +
    # DMRs
    geom_rect(data = dmrPlot %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = (as.integer(chr) - 0.1 - percMax), xmax = as.integer(chr) - 0.1, ymin = pos, ymax = pos+1),
              fill = "#ae2012", size = 0.25) +
    # centromeres
    geom_point(data = centPos %>% filter(chr %in% paste0("chr", 1:12)), aes(x = chr, y = centre), 
               size = 5, colour = "black") +
    # legend bars
    geom_rect(data = adL, 
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              size = 0.25) +
    # legend text
    geom_text(data = adW, 
              aes(x = x, y = y, label = label))+
    scale_fill_manual(values = c("#005f73","#ae2012")) +
    ylim(0, 250) +
    theme_classic() +
    theme(text = element_text(size=15),axis.line=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")+
    labs(x=NULL,y="Chromosome Size (Mb)")
  
  # chromosomes 13 - 22 + X
  
  # very annoying but you have to filter all the dataframes or else the factor levels won't match the integer value
  chromSizeFilt <- chromSize %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  chromSizeFilt$chr <- factor(chromSizeFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  centPosFilt <- centPos %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  centPosFilt$chr <- factor(centPosFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  asePlotFilt <- asePlot %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  asePlotFilt$chr <- factor(asePlotFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  dmrPlotFilt <- dmrPlot %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  dmrPlotFilt$chr <- factor(dmrPlotFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  
  # chromosomes 13-22+X
  p2 <- ggplot() +
    # chromosome bars
    geom_segment(data = chromSizeFilt, aes(x = chr, xend = chr, y = 0, yend = size), 
                 lineend = "round", color = "lightgrey", size = 5) +
    # ASE genes
    geom_rect(data = asePlotFilt, 
              aes(xmin = as.integer(chr) + 0.1, xmax = (as.integer(chr) + 0.1 + percMax), ymin = pos, ymax = pos+1),
              fill = "#005f73", size = 0.25) +
    # DMRs
    geom_rect(data = dmrPlotFilt, 
              aes(xmin = (as.integer(chr) - 0.1 - percMax), xmax = as.integer(chr) - 0.1, ymin = pos, ymax = pos+1),
              fill = "#ae2012", size = 0.25) +
    # centromeres
    geom_point(data = centPosFilt, aes(x = chr, y = centre), 
               size = 5, colour = "black") +
    ylim(0, 250) +
    theme_classic() +
    theme(text = element_text(size=15),axis.line=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")+
    labs(x=NULL,y=NULL)
  
} else{ ## NO DMRs OR CNVs
  # legend
  adL <- data.frame(xmin = c(9.7, 9.7, 10.3), xmax = c(10.35,9.75,10.35), ymin = c(240,238,238), ymax = c(243,243,243),
                    fill = c("ase","ase","ase"))
  
  adW <- data.frame(x = c(11.5,10), y = c(240,250),
                    label = c("ASE Gene Density", as.character(c(maxASE))))
  
  # plot
  # chromosomes 1 - 12
  p1 <- ggplot() +
    # chromosome bars
    geom_segment(data = chromSize %>% filter(chr %in% paste0("chr", 1:12)), aes(x = chr, xend = chr, y = 0, yend = size), 
                 lineend = "round", color = "lightgrey", size = 5) +
    # ASE genes
    geom_rect(data = asePlot %>% filter(chr %in% paste0("chr", 1:12)), 
              aes(xmin = as.integer(chr) + 0.1, xmax = (as.integer(chr) + 0.1 + percMax), ymin = pos, ymax = pos+1),
              fill = "#005f73", size = 0.25) +
    # centromeres
    geom_point(data = centPos %>% filter(chr %in% paste0("chr", 1:12)), aes(x = chr, y = centre), 
               size = 5, colour = "black") +
    # legend bars
    geom_rect(data = adL, 
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              size = 0.25) +
    # legend text
    geom_text(data = adW, 
              aes(x = x, y = y, label = label))+
    scale_fill_manual(values = c("#005f73")) +
    ylim(0, 250) +
    theme_classic() +
    theme(text = element_text(size=15),axis.line=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")+
    labs(x=NULL,y="Chromosome Size (Mb)")
  
  # chromosomes 13 - 22 + X
  
  # very annoying but you have to filter all the dataframes or else the factor levels won't match the integer value
  chromSizeFilt <- chromSize %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  chromSizeFilt$chr <- factor(chromSizeFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  centPosFilt <- centPos %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  centPosFilt$chr <- factor(centPosFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  asePlotFilt <- asePlot %>% filter(chr %in% c(paste0("chr", 13:22), "chrX"))
  asePlotFilt$chr <- factor(asePlotFilt$chr,levels=c(paste0("chr", 13:22), "chrX"))
  
  # chromosomes 13-22+X
  p2 <- ggplot() +
    # chromosome bars
    geom_segment(data = chromSizeFilt, aes(x = chr, xend = chr, y = 0, yend = size), 
                 lineend = "round", color = "lightgrey", size = 5) +
    # ASE genes
    geom_rect(data = asePlotFilt, 
              aes(xmin = as.integer(chr) + 0.1, xmax = (as.integer(chr) + 0.1 + percMax), ymin = pos, ymax = pos+1),
              fill = "#005f73", size = 0.25) +
    # centromeres
    geom_point(data = centPosFilt, aes(x = chr, y = centre), 
               size = 5, colour = "black") +
    ylim(0, 250) +
    theme_classic() +
    theme(text = element_text(size=15),axis.line=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = "none")+
    labs(x=NULL,y=NULL)
}


## ---------------------------------------------------------------------------
## PLOT 
## ---------------------------------------------------------------------------

# put them together
plot <- plot_grid(p1, p2, align = "v", axis = "l", nrow = 2)

# save plot
ggsave(plot, filename = paste0(opt$out,".pdf"), width = 10, height = 7, units = "in")

