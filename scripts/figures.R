
## ---------------------------------------------------------------------------
## Allelic Imbalance in Expression using MBASED, part 2
## Vanessa Porter, Mar. 2022
## ---------------------------------------------------------------------------

suppressMessages(library(optparse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tibble))
suppressMessages(library(chromPlot))
suppressMessages(library(networkD3))
suppressMessages(library(htmlwidgets))
suppressMessages(library(reshape2))
suppressMessages(library(ggrepel))
suppressMessages(library(ggsci))

## ---------------------------------------------------------------------------
## FIGURES
## ---------------------------------------------------------------------------


# Make help options
option_list = list(
  make_option(c("-b", "--mbased"), type="character", default=NULL,
              help="mbased dataframe output", metavar="character"),
  make_option(c("-r", "--rpkm"), type="character", default=NULL,
              help="RPKM matrix", metavar="character"),
  make_option(c("-g", "--gene"), type="character", default = NULL,
              help="Ensembl gene annotation", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default = NULL,
              help="Sample name", metavar="character"),
  make_option(c("-m", "--min"), type="numeric", default = 1,
              help="Minimum RPKM value", metavar="numeric"),
  make_option(c("-o", "--outdir"), type="character", default = NULL,
              help="Output directory", metavar="character")
)

### 
### SET UP THE DATAFRAME
###

# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
df <- read.delim(opt$mbased, header = T, stringsAsFactors = F)
rpkm <- read.delim(opt$rpkm, header = T, stringsAsFactors = F)
all_genes <- read.delim(opt$gene, header = F, stringsAsFactors = F)
out <- opt$outdir
sample <- opt$sample
min <- opt$min

# fix sample name
sample <- ifelse(length(grep("-", sample)) == 0, sample, gsub("-", ".", sample))

# select sample
rpkm_sample <- rpkm[,c("gene", sample)] 
colnames(rpkm_sample) <- c("gene", "expr")

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

print("Beginning figures ...")

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
#### SANKEY PLOT
####

# set filters on the RPKM matrix
rpkm_sample$gene_biotype <- all_genes$V7[match(rpkm_sample$gene, all_genes$V4)]
rpkm_sample_filt1 <- rpkm_sample[rpkm_sample$gene_biotype %in% c("lincRNA", "miRNA", "protein_coding"),]
rpkm_sample_filt2 <- rpkm_sample_filt1[rpkm_sample_filt1$expr > min,] 

# get the input values for the plot
a <- nrow(rpkm_sample_filt1)
b <- c(nrow(rpkm_sample_filt2),nrow(rpkm_sample_filt1[rpkm_sample_filt1$expr <= min,] ))
c <- c(nrow(df),nrow(rpkm_sample_filt2[!rpkm_sample_filt2$gene %in% df$gene,]))
d <- c(nrow(df[df$padj < 0.05 & df$majorAlleleFrequency > 0.75,]),
       sum(nrow(df[df$padj >= 0.05 & df$majorAlleleFrequency <= 0.75,]),
           nrow(df[df$padj >= 0.05 & df$majorAlleleFrequency > 0.75,]),
           nrow(df[df$padj < 0.05 & df$majorAlleleFrequency <= 0.75,])))

# create a connection data frame
links <- data.frame(
  source=c(rep(paste0("All Genes (n=", a, ")"), 2), rep(paste0("Expressed (n=", b[1], ")"), 2), rep(paste0("Phased Genes (n=", c[1], ")"), 2)),
  target=c(paste0("Expressed (n=", b[1], ")"), paste0("Not Expressed (n=", b[2], ")"), 
           paste0("Phased Genes (n=", c[1], ")"), paste0("Unphased Genes (n=", c[2], ")"), 
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

saveWidget(sankey, file=paste0(out, "/sankeyPlot.html"), selfcontained = F)

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


print("Figures completed")
