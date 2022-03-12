
# Allele specific expression
This workflow outputs the allele-specific expression using short-read RNA-seq and DNA-seq bams. The phasing information of the variants from tools such as WhatsApp can be provided to increase the performance of the tool (to-be an optional input but I don't know how to write that in snakemake - it is optional in the R script though)

# Install
Clone the repository:

```
git clone https://svn.bcgsc.ca/bitbucket/scm/marra/vporter-allelespecificexpression.git
```

# Running samples
## **Edit the config files**

parameters.yaml: <br />
- Choose a genome to use (hg38/hg19/hg38_no_alt_TCGA_HTMCP_HPVs)
- Specify the path to the RPKM matrix that contains the sample (can include multiple samples but sample id must be in the column names)

samples.yaml: <br />
- Make sure to name the sample(s) the same identifier as the matching RPKM matrix column name
- Add the paths to the DNA-seq and RNA-seq Illumina bam files and the phased VCF file 

## **Run snakemake**
Right now conda environments are used, but it would be good to eventually change this to docker containers. You can choose the number of max threads to use with `-c`. This is the command to run:

```
snakemake --use-conda -c 30
```

# Outputs
## MBASED-related output files
All of the outputs will be found in the `output/{sample}` directory, but the main outputs of interest will be in `output/{sample}/mbased`. The outputs include:
1. The tabular results of the output `MBASED_expr_gene_results.txt`
2. The rds object of the MBASED raw output `MBASEDresults.rds`
3. The rds object of the MBASED gene output simplified dataframe `MBASED_expr_gene_results.rds`
4. Four descriptive figures of the results: `aseGenesBar.pdf`, `aseGenesDot.pdf`, `chromPlot.pdf`, `sankeyPlot.html` 

## Output dataframe description file 
The `MBASED_expr_gene_results.txt` output will have the main results included. Here is a description of the columns:
 <br />

| Column               | Description                                                           | 
| :---                 |    :----:                                                             |  
| gene                 | HGNC gene symbol                                                      | 
| gene_biotype         | The type of gene (protein_coding currently selected)                  | 
| gene_band            | chromosome band                                                       | 
| RPKM                 | Expression level                                                      | 
| allele1IsMajor       | T/F is allele 1 is the major allele (allele 1 = HP1)                  | 
| majorAlleleFrequency | Major allele frequency                                                | 
| pValueASE            | p-value output by mBASED                                              | 
| pValueHeterogeneity  | The heterogeneity of the pvalue between the multiple SNVs in the gene | 
| padj                 | Benjamini-Hochberg adjusted pvalue                                    | 
| significance         | Factorize p-value cut off of 0.05                                     | 
| MAF                  | Factorize MAF cut off of 0.75                                         | 

## Figure Outputs 

### aseGenesDot
The `aseGenesDot.pdf` figure shows the MAF and padj for each phased gene. The red dots (MAF > 0.75 & padj < 0.05) are the ASE genes. 

![Allele Specific Expression Bar Graph](example_figures/aseGenesDot.png "ASE Dot Plot")

### aseGenesBar
The `aseGenesBar.pdf` figure shows the number of genes that fall into the statistical cut offs across each chromosome. The red variable (MAF > 0.75 & padj < 0.05) are the ASE genes. The remaining genes are coloured by where they fall in the MAF and padj cut-offs. 

![Allele Specific Expression Bar Graph](example_figures/aseGenesBar.png "ASE Bar Graph")

### chromPlot
The `chromPlot.pdf` figure shows the distribution of ASE genes by their genomic location. The ASE genes are in red facing left, while the biallelic genes (all genes that don't make the padj and MAF cut offs) are in grey facing right. 

![Allele Specific Expression Bar Graph](example_figures/chromPlot.png "ASE Chromosome Plot")

### sankeyPlot
The `sankeyPlot.html` shows how many of genes made it through each filtering step. All genes included only the genes labeled as protein coding, lincRNA, or miRNA.  They were filtered by (1) expression level, (2) ablity to be phased, and (3) whether or not they have ASE.

![Allele Specific Expression Bar Graph](example_figures/sankeyPlot.png "ASE Sankey Plot")