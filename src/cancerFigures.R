suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library("ggsci"))


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
summary_table <- read.delim(summary_path, header = T, stringsAsFactors = F, sep = "\t")


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
    xlab("Copy number variant state") + 
    theme_bw() + 
    scale_color_npg() + 
    scale_fill_npg() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  
  
  ggsave(filename = paste0(out, "/cnvBar.pdf"), plot = cnvBar, units = "in")
  
  
  expectMAF <- summary_table %>%
    dplyr::filter(!is.na(cnv_state)) %>%
    dplyr::mutate(cnvRatioDiff = majorAlleleFrequency - expectedMAF) %>%
    ggplot(aes(cnvRatioDiff, padj)) + 
    geom_point(aes(colour = cnv_state)) +
    geom_vline(xintercept = 0.10) +
    geom_hline(yintercept = 0.05) +
    xlab("MAF - Expectd MAF from CNV") + 
    ylab("Adjusted pvalue") + 
    ggtitle(paste0("Copy Number Variants adjusted for expected MAF"), 
            subtitle = sample) +
    theme_bw() + 
    scale_fill_npg() +
    scale_color_npg() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
  
  ggsave(filename = paste0(out, "/expectMAF.pdf"), plot = expectMAF, units = "in")
  
  #### TEST CORRELATION PLOT ####
  correlationPlot <- summary_table %>%
    dplyr::filter(cnv_state == "imbalance") %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::filter(!is.na(cnv_state)) %>%
    dplyr::mutate(rawExpectedMAF = pmax(cnv.A, cnv.B)/(cnv.A + cnv.B)) %>%
    ggplot(aes(majorAlleleFrequency, rawExpectedMAF)) +
    geom_point(aes(color = aseResults), alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE) +
    ylab("CNV expected MAF") + 
    xlab("MBASED Major Allele Frequency") + 
    ggtitle("Correlation between CNV expected MAF and \nMBASED MAF for CNV Imbalanced genes",
            subtitle = sample) + 
    theme_bw() + 
    scale_fill_npg() +
    scale_color_npg() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) 
    ggsave(filename = paste0(out, "/cnvMAFcorrlation.pdf"), plot = correlationPlot, units = "in")
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
      methyl_state < 0 ~ "allele2Methyl",
      TRUE ~ "allele1Methyl"
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
    ggtitle("Allelic Methylation and Expression Contigency Table", subtitle = sample) +
    theme_bw() + 
    theme(legend.position = "none", 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(), 
          axis.title.x=element_blank(), 
          axis.title.y=element_blank())  
  
  ggsave(filename = paste0(out, "/dmrContingency.pdf"), plot = dmr_contigency, units = "in")
  
  dmr %>%
    dplyr::mutate(sample = sample) %>%
    write.table(file = paste0(out, "/tables/dmrContingency.tsv"),
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
    write.table(file = paste0(out, "/tables/stopVarContingency.tsv"),
                quote = F, sep = "\t", row.names = T, col.names = T)
  
  stopContingency <- ggplot(data =  stop, mapping = aes(x = stop, y = expression)) +
      geom_tile(aes(fill = n), colour = "white") +
      geom_text(aes(label = sprintf("%1.0f", n)), vjust = 1, colour = "white") +
      scale_fill_gradient(low = "blue", high = "red") +
      theme_bw() +
      ggtitle("ASE vs Stop gain/loss contingency table", subtitle = sample) + 
      theme(legend.position = "none", 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.border = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(), 
            axis.title.x=element_blank(), 
            axis.title.y=element_blank()) 
  ggsave(filename = paste0(out, "/stopVarContingency.pdf"), plot = stopContingency, units = "in")
  
  stopVarTable <- summary_table %>%
    dplyr::mutate(stopVar = !is.na(stop_variant_allele)) %>%
    dplyr::select(aseResults, stopVar) %>%
    group_by(aseResults, stopVar) %>%
    summarize(n=n()) %>%
    dplyr::mutate(percent = n/sum(n)) %>%
    dplyr::mutate(sample = sample)
  
  write.table(stopVarTable, file = paste0(out, "/tables/stopVarTable.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  
  stopBar <-  stopVarTable %>%
    dplyr::filter(stopVar) %>%
    ggplot(aes(aseResults, percent)) +
    geom_bar(aes(fill = aseResults ), stat = "identity") +
    ylab("Proportion of genes") +
    coord_flip() +
    ggtitle("Proportion of genes with stop gain or loss variant", 
            subtitle = sample) + 
    theme_bw() + 
    scale_fill_npg() +
    scale_color_npg() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          legend.position = "none") 
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
  
  write.table(snvTable, file = paste0(out, "/tables/snvTable.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  
  snv <- snvTable %>%
    dplyr::filter(somaticSNV) %>%
    ggplot(aes(aseResults, percent)) +
    geom_bar(aes(fill = aseResults), stat = "identity") +
    ylab("Proportion of genes") +
    coord_flip() +
    ggtitle("Proportion of genes with somatic SNV mutation", 
            subtitle = sample) + 
    theme_bw() + 
    scale_fill_npg() +
    scale_color_npg() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          legend.position = "none") 
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
  
  write.table(indelTable, file = paste0(out, "/tables/indelTable.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  
  indel <- indelTable %>%
    dplyr::filter(somaticIndel) %>%
    ggplot(aes(aseResults, percent)) +
    geom_bar(aes(fill = aseResults), stat = "identity") +
    ylab("Proportion of genes") +
    coord_flip() +
    ggtitle("Proportion of genes with somatic indel mutation", 
            subtitle = sample) + 
    theme_bw() + 
    scale_fill_npg() +
    scale_color_npg() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          legend.position = "none") 
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
  write.table(tfbsTable, file = paste0(out, "/tables/tfbsTable.tsv"),
              quote = F, sep = "\t", row.names = F, col.names = T)
  
  tfbs <- tfbsTable %>%
    dplyr::filter(tfbsVar) %>%
    ggplot(aes(aseResults, percent)) +
    geom_bar(aes(fill = aseResults), stat = "identity")+
    ylab("Proportion of genes") +
    coord_flip() +
    ggtitle("Proportion of genes with mutation in \ntranscription factor binding site", 
            subtitle = sample) + 
    theme_bw() + 
    scale_fill_npg() +
    scale_color_npg() + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(), 
          legend.position = "none") 
  ggsave(filename = paste0(out, "/tfbsBar.pdf"), plot = tfbs, units = "in")
}


#######################
# Cancer Figures
#######################
cancer_bar <- summary_table %>%
  dplyr::filter(cancer_gene) %>%
  ggplot(aes(aseResults)) + 
  geom_bar(aes(fill = aseResults)) +
  theme_bw() +
  scale_fill_npg() +
  scale_color_npg() +
  ggtitle("MBASED result for cancer genes", subtitle=sample) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = "none") 
ggsave(filename = paste0(out, "/cancerBar.pdf"), plot = cancer_bar, units = "in")

#######################
# Mechanism function
#######################

causeFunction <- function(dataset) {
  summaryTableExplain <- dataset %>%
    dplyr::filter(aseResults == "ASE") 
  ASEunknown <- summaryTableExplain
  cause_final <- tibble()
  
  if ("cnv_state" %in% column) {
    
    loh <- summaryTableExplain %>%
      dplyr::filter(cnv_state == "LOH") %>%
      nrow()
    cause_final <- bind_rows(cause_final, tibble(cause = "LOH", num = loh))
    
    ASEunknown <- ASEunknown %>%
      dplyr::filter(cnv_state != "LOH")
    
    summaryTableExplain$cnvRatioDiff = summaryTableExplain$majorAlleleFrequency - summaryTableExplain$expectedMAF
    ASEunknown$cnvRatioDiff = ASEunknown$majorAlleleFrequency - ASEunknown$expectedMAF
    
    cnv <- summaryTableExplain %>%
      dplyr::filter(cnv_state != "LOH") %>%
      dplyr::filter(cnvRatioDiff< 0.10) %>%
      nrow()
    cause_final <- bind_rows(cause_final, tibble(cause = "CNV imbalance", num = cnv))
    ASEunknown <- ASEunknown %>%
      dplyr::filter(cnvRatioDiff <= 0.10)
  }
  
  if ("somaticSNV" %in% column) {
    somaticSNV <- summaryTableExplain %>%
      dplyr::filter(somaticSNV) %>%
      nrow()
    cause_final <- bind_rows(cause_final, tibble(cause = "Somatic SNV", num = somaticSNV))
    ASEunknown <- ASEunknown %>%
      dplyr::filter(!somaticSNV)
  }
  if ("somaticIndel" %in% column) {
    somaticIndel <- summaryTableExplain %>%
      dplyr::filter(somaticIndel) %>%
      nrow()
    cause_final <- bind_rows(cause_final, tibble(cause = "Somatic Indel", num = somaticIndel))
    ASEunknown <- ASEunknown %>%
      dplyr::filter(!somaticIndel)
  }
  
  if ("methyl_state" %in% column) {
    methyl <- summaryTableExplain %>%
      dplyr::filter((methyl_state < 0 & allele1IsMajor) | (methyl_state > 0 & !allele1IsMajor)) %>%
      nrow()
    cause_final <- bind_rows(cause_final, tibble(cause = "Allelic methylation", num = methyl))
    ASEunknown <- ASEunknown %>%
      dplyr::filter(is.na(methyl_state) | !((methyl_state < 0 & allele1IsMajor) | (methyl_state > 0 & !allele1IsMajor)))
  }
  
  if ("stop_variant_allele" %in% column) {
    stop <- summaryTableExplain %>% 
      dplyr::filter((allele1IsMajor & stop_variant_allele == 2) | (!allele1IsMajor & stop_variant_allele == 1)) %>%
      nrow()
    cause_final <- bind_rows(cause_final, tibble(cause = "Stop variant", num = stop))
    ASEunknown <- ASEunknown %>%
      dplyr::filter(is.na(stop_variant_allele) | !((allele1IsMajor & stop_variant_allele == 2) | (!allele1IsMajor & stop_variant_allele == 1)))
  }
  if ("tf_allele" %in% column) {
    tfbs <- summaryTableExplain %>%
      dplyr::filter(!is.na(transcriptionFactor)) %>%
      nrow()
    cause_final <- bind_rows(cause_final, tibble(cause = "TFBS variant", num = tfbs))
    ASEunknown <- ASEunknown %>%
      dplyr::filter(is.na(transcriptionFactor))
  }
  
  unknown <- nrow(ASEunknown)
  cause_final <- bind_rows(cause_final, tibble(cause = "Unknown", num = unknown))

  
   return(cause_final)
}

################
# ASE mechanism
################
aseMechanism <- causeFunction(summary_table)

aseCause <- aseMechanism %>%
  ggplot(aes(reorder(cause, -num), num)) + 
  geom_bar(aes(fill = cause), stat = "identity") +
  theme_bw() + 
  scale_fill_npg() +
  scale_color_npg() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  xlab("ASE mechanism") + 
  ylab("Gene frequency") + 
  ggtitle("Frequency of all ASE genes for each genetic mechanism", 
          subtitle = sample)

ggsave(filename = paste0(out, "/aseCause.pdf"), plot = aseCause, units = "in")

write.table(aseMechanism, file = paste0(out, "/tables/aseCause.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)

################################
# Cancer ASE mechanism
################################
cancerMechanism <- summary_table %>%
  dplyr::filter(cancer_gene) %>%
  causeFunction()


cancerCause <- cancerMechanism %>%
  ggplot(aes(reorder(cause, -num), num)) + 
  geom_bar(aes(fill = cause), stat = "identity") +
  theme_bw() + 
  scale_fill_npg() +
  scale_color_npg() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(), 
        legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  xlab("ASE mechanism") + 
  ylab("Gene frequency") + 
  ggtitle("Frequency of cancer ASE genes for each genetic mechanism", 
          subtitle = sample)
ggsave(filename = paste0(out, "/cancerCause.pdf"), plot = cancerCause, units = "in")


write.table(cancerMechanism, file = paste0(out, "/tables/cancerCause.tsv"),
            quote = F, sep = "\t", row.names = F, col.names = T)
