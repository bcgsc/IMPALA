suppressMessages(library(tidyverse))
suppressMessages(library(optparse))


# Make help options
option_list = list(
  make_option(c("-s", "--summary"), type="character", default=NULL,
              help="Summary Table file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default = NULL,
              help="Output directory name", metavar="character"),
  make_option(c("-a", "--sample"), type="character", default = NULL,
              help="Sample name", metavar="character")
)


# load in options 
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
out <- opt$outdir
sample <- opt$sample

summary_path <- opt$summary
summary_table <- read_delim(summary_path, delim = "\t", show_col_types = F)

column <- colnames(summary_table)

############
#CNV Figures
############
if ("cnv_state" %in% column) {
  print("Creating CNV figure...")
  cnvBar <- summary_table %>%
    dplyr::filter(!is.na(cnv_state)) %>%
    ggplot(aes(cnv_state)) + 
    geom_bar(aes(fill = aseResults), position = "dodge") +
    ggtitle(paste0("Gene frequency in each copy number variant state"), 
            subtitle = sample) +
    ylab("Gene Count") + 
    xlab("Copy number variant state") 
  
  ggsave(filename = paste0(out, "/cnvBar.pdf"), plot = cnvBar, units = "in")
  
  
  expectMAF <- summary_table %>%
    dplyr::filter(!is.na(cnv_state)) %>%
    dplyr::mutate(cnvRatioDiff = majorAlleleFrequency - expectedMAF) %>%
    ggplot(aes(cnvRatioDiff, padj)) + 
    geom_point(aes(colour = cnv_state)) +
    geom_vline(xintercept = 0.10) +
    geom_hline(yintercept = 0.05) +
    xlab("MAF - Expectd MAF from CNV") + 
    ylab("Adjusted pvalue")
  
  ggsave(filename = paste0(out, "/expectMAF.pdf"), plot = expectMAF, units = "in")
  
}
 

##############
# DMR Figures
##############

if ("methyl_state" %in% column) {
  print("Creating DMR figures...")
  
  dmr <- summary_table %>%
    dplyr::filter(!is.na(methyl_state)) %>%
    dplyr::filter(aseResults == "ASE") %>%
    dplyr::mutate(methylation = case_when(
      methyl_state < 0 ~ "allele1Methyl",
      TRUE ~ "allele2Methyl"
    )) %>%
    dplyr::mutate(expression = case_when(
      allele1IsMajor ~ "alelle1Expression",
      TRUE ~ "alelle2Expression"
    )) %>%
    dplyr::select(methylation, expression) %>%
    group_by(methylation, expression) %>%
    summarize(n=n()) 
    
  
  dmr_contigency <- ggplot(data =  dmr, mapping = aes(x = methylation, y = expression)) +
    geom_tile(aes(fill = n), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, colour = "white") +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_bw() + theme(legend.position = "none")
  ggsave(filename = paste0(out, "/dmrContingency.pdf"), plot = dmr_contigency, units = "in")
  
  dmr %>%
    dplyr::mutate(sample = sample) %>%
    write.table(file = paste0(out, "/table/dmrContingency.tsv"),
                quote = F, sep = "\t", row.names = T, col.names = T)
}




##############
# Stop Var Figures
##############
if ("stop_variant_allele" %in% column) {
  print("Creating Stop gain/loss figures...")
  
  stop <- summary_table %>%
    dplyr::select(gene, allele1IsMajor, aseResults, stop_variant_allele) %>% 
    dplyr::filter(!is.na(stop_variant_allele)) %>%
    dplyr::filter(aseResults == "ASE") %>%
    dplyr::mutate(expression = case_when(
      allele1IsMajor ~ "alelle1Expression",
      TRUE ~ "alelle2Expression"
    )) %>%
    dplyr::mutate(stop = case_when(
      stop_variant_allele == 1 ~ "alelle1StopVar",
      TRUE ~ "alelle2StopVar"
    )) %>%
    group_by(expression, stop) %>%
    summarize(n=n()) 
  
  stop %>%
    dplyr::mutate(sample = sample) %>%
    write.table(file = paste0(out, "/table/stopVarContingency.tsv"),
                quote = F, sep = "\t", row.names = T, col.names = T)
  
  stopContingency <- ggplot(data =  stop, mapping = aes(x = stop, y = expression)) +
      geom_tile(aes(fill = n), colour = "white") +
      geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, colour = "white") +
      scale_fill_gradient(low = "blue", high = "red") +
      theme_bw() + theme(legend.position = "none") +
      ggtitle("ASE vs Stop gain/loss contingency table", subtitle = sample)
  ggsave(filename = paste0(out, "/stopVarContingency.pdf"), plot = stopContingency, units = "in")
  
  stopVarTable <- summary_table %>%
    dplyr::mutate(stopVar = !is.na(stop_variant_allele)) %>%
    dplyr::select(aseResults, stopVar) %>%
    group_by(aseResults, stopVar) %>%
    summarize(n=n()) %>%
    dplyr::mutate(percent = n/sum(n)) %>%
    dplyr::mutate(sample = sample)
  
  write.table(stopVarTable, file = paste0(out, "/table/stopVarTable.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  
  stopBar <-  stopVarTable %>%
    dplyr::filter(stopVar) %>%
    ggplot(aes(aseResults, percent)) +
    geom_bar(stat = "identity") +
    ylab("Proportion of genes") +
    coord_flip() +
    ggtitle("Proportion of genes with stop gain or loss variant", 
            subtitle = sample)
  ggsave(filename = paste0(out, "/stopVarBar.pdf"), plot = stopBar, units = "in")
}
  
  


##############
# Somatic SNV
##############
if ("somaticSNV" %in% column) {
  print("Creating somatic SNV figures...")
  
  snvTable <- summary_table %>%
    dplyr::select(gene, aseResults, somaticSNV) %>%
    group_by(aseResults, somaticSNV) %>%
    summarize(n=n()) %>%
    dplyr::mutate(percent = n/sum(n)) %>%
    dplyr::mutate(sample = sample)
  
  write.table(snvTable, file = paste0(out, "/table/snvTable.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  
  snv <- snvTable %>%
    dplyr::filter(somaticSNV) %>%
    ggplot(aes(aseResults, percent)) +
    geom_bar(stat = "identity") +
    ylab("Proportion of genes") +
    coord_flip() +
    ggtitle("Proportion of genes with somatic SNV mutation", 
            subtitle = sample)
  ggsave(filename = paste0(out, "/somaticSNVbar.pdf"), plot = snv, units = "in")
}



##############
# Somatic Indel
##############
if ("somaticIndel" %in% column) {
  print("Creating somatic indel figures...")
  
  indelTable <- summary_table %>%
    dplyr::select(gene, aseResults, somaticIndel) %>%
    group_by(aseResults, somaticIndel) %>%
    summarize(n=n()) %>%
    dplyr::mutate(percent = n/sum(n)) %>%
    dplyr::mutate(sample = sample)
  
  write.table(indelTable, file = paste0(out, "/table/indelTable.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  
  indel <- indelTable %>%
    dplyr::filter(somaticIndel) %>%
    ggplot(aes(aseResults, percent)) +
    geom_bar(stat = "identity") +
    ylab("Proportion of genes") +
    coord_flip() +
    ggtitle("Proportion of genes with somatic indel mutation", 
            subtitle = sample)
  ggsave(filename = paste0(out, "/somaticIndelBar.pdf"), plot = indel, units = "in")
}


##############
# TFBS Figures
##############
if ("tf_allele" %in% column) {
  print("Creating TFBS figures...")
  
  tfbsTable <- summary_table %>%
    dplyr::select(gene, aseResults, transcriptionFactor) %>%
    dplyr::mutate(tfbsVar = !is.na(transcriptionFactor)) %>%
    group_by(aseResults, tfbsVar) %>%
    summarize(n=n()) %>%
    dplyr::mutate(percent = n/sum(n)) %>%
    dplyr::mutate(sample = sample)
  write.table(tfbsTable, file = paste0(out, "/table/tfbsTable.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  
  tfbs <- tfbsTable %>%
    dplyr::filter(tfbsVar) %>%
    ggplot(aes(aseResults, percent)) +
    geom_bar(stat = "identity")+
    ylab("Proportion of genes") +
    coord_flip() +
    ggtitle("Proportion of genes with mutation in transcription factor binding site", 
            subtitle = sample)
  ggsave(filename = paste0(out, "/tfbsBar.pdf"), plot = tfbs, units = "in")
}


#######################
# Explanation Pie Chart
#######################

summaryTableExplain <- data.frame(summary_table) %>%
  dplyr::filter(aseResults == "ASE") %>%
  dplyr::mutate(explanation = NA)

if ("cnv_state" %in% column) {
  summaryTableExplain$cnvRatioDiff = summaryTableExplain$majorAlleleFrequency - summaryTableExplain$expectedMAF
  summaryTableExplain$explanation[is.na(summaryTableExplain$explanation) & 
                                summaryTableExplain$cnv_state=="LOH"] <- "LOH"
  summaryTableExplain$explanation[is.na(summaryTableExplain$explanation) & 
                                summaryTableExplain$cnvRatioDiff< 0.10] <- "CNV"
}

if ("somaticSNV" %in% column) {
  summaryTableExplain$explanation[is.na(summaryTableExplain$explanation) & 
                                summaryTableExplain$somaticSNV] <- "Somatic Mutation"
}
if ("somaticIndel" %in% column) {
  summaryTableExplain$explanation[is.na(summaryTableExplain$explanation) & 
                                summaryTableExplain$somaticIndel] <- "Somatic Mutation"
}

if ("methyl_state" %in% column) {
  summaryTableExplain$explanation[is.na(summaryTableExplain$explanation) & 
                                summaryTableExplain$allele1IsMajor & 
                                !is.na(summaryTableExplain$methyl_state) & 
                                summaryTableExplain$methyl_state < 0] <- "Alleleic Methylation"
  summaryTableExplain$explanation[is.na(summaryTableExplain$explanation) & 
                                !summaryTableExplain$allele1IsMajor & 
                                !is.na(summaryTableExplain$methyl_state) & 
                                summaryTableExplain$methyl_state > 0] <- "Alleleic Methylation"
}

if ("stop_variant_allele" %in% column) {
  summaryTableExplain$explanation[is.na(summaryTableExplain$explanation) & 
                                summaryTableExplain$allele1IsMajor & 
                                summaryTableExplain$stop_variant_allele == 2] <- "Stop gain/loss mutation"
  summaryTableExplain$explanation[is.na(summaryTableExplain$explanation) & 
                                !summaryTableExplain$allele1IsMajor & 
                                summaryTableExplain$stop_variant_allele == 1] <- "Stop gain/loss mutation"
}
if ("tf_allele" %in% column) {
  summaryTableExplain$explanation[is.na(summaryTableExplain$explanation) & 
                                !is.na(summaryTableExplain$transcriptionFactor)] <- "TFBS mutation"
}


summaryTableExplain$explanation[is.na(summaryTableExplain$explanation)] <- "Unknown"

summaryTableExplain %>%
  dplyr::select(gene, cancer_gene, explanation) %>%
  dplyr::mutate(sample = sample) %>%
  write.table(file = paste0(out, "/table/aseCause.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)


asePie <- summaryTableExplain %>%
  group_by(explanation) %>%
  summarize(n=n()) %>%
  ggplot(aes(x="", y=n, fill=explanation)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  ggtitle("Causes of allelic expression", subtitle = sample) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank()) + 
  labs(fill='ASE Causes') 

ggsave(filename = paste0(out, "/aseCause.pdf"), plot = asePie, units = "in")

#######################
# Cancer Figures
#######################
cancer_bar <- summary_table %>%
  dplyr::filter(cancer_gene) %>%
  ggplot(aes(aseResults)) + 
  geom_bar() +
  ggtitle("MBASED result for cancer genes", subtitle=sample)
ggsave(filename = paste0(out, "/cancerBar.pdf"), plot = cancer_bar, units = "in")

#######################
# Cancer Explain Figures
#######################
cancerCause <- summaryTableExplain %>%
  dplyr::filter(cancer_gene) %>%
  group_by(explanation) %>%
  summarize(n=n()) %>%
  ggplot(aes(x="", y=n, fill=explanation)) +
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) +
  ggtitle("Cause of allelic expression in cancer genes", subtitle = sample) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks = element_blank())

ggsave(filename = paste0(out, "/cancerCause.pdf"), plot = cancerCause, units = "in")
