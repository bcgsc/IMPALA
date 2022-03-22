####
#### FILTER THE ENSEMBL BIOMART FOR ASE
####
suppressMessages(library(biomaRt))

# select mart
grch38 <- useEnsembl(biomart = 'genes', 
           dataset = 'hsapiens_gene_ensembl',
           version = 100)

# retrieve gene info from biomaRt
all_genes <- getBM(attributes = c("chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol", "band", "gene_biotype"),mart = grch38)

# filter and adjust chr name
all_genes <- all_genes_nofilt[all_genes_nofilt$gene_biotype %in% c("lincRNA", "miRNA", "protein_coding"),]
all_genes <- all_genes[all_genes$chromosome_name %in% c(as.character(1:22), "X", "Y"),]
all_genes <- all_genes[all_genes$hgnc_symbol != "",]
all_genes$locus <- paste0(all_genes$chromosome_name, all_genes$band)
all_genes$chromosome_name <- paste0("chr", all_genes$chromosome_name)

# save
write.table(all_genes, file="/projects/vporter_prj/git-repo/ASE/annotation/biomart_ensembl100_GRCh38.bed", quote = F, sep = "\t", col.names = F, row.names = F)

