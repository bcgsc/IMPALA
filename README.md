
## Allele specific expression
This workflow outputs the allele-specific expression using short-read RNA-seq and DNA-seq bams. The phasing information of the variants from tools such as WhatsApp can be provided to increase the performance of the tool (to-be an optional input but I don't know how to write that in snakemake - it is optional in the R script though)

### Install
Clone the repository:

```
git clone https://svn.bcgsc.ca/bitbucket/scm/marra/vporter-allelespecificexpression.git
```

### Running samples
1. **Edit the config files** <br />
parameters.yaml: <br />
- Choose a genome to use (hg38/hg19/hg38_no_alt_TCGA_HTMCP_HPVs)
- Specify the path to the RPKM matrix that contains the sample (can include multiple samples but sample id must be in the column names)

samples.yaml: <br />
- Make sure to name the sample(s) the same identifier as the matching RPKM matrix column name
- Add the paths to the DNA-seq and RNA-seq Illumina bam files and the phased VCF file 

2. **Run snakemake** <br />
Right now conda environments are used, but it would be good to eventually change this to docker containers. You can choose the number of max threads to use with `-c`. This is the command to run:

```
snakemake --use-conda -c 30
```

### Outputs
All of the outputs will be found in the `output/{sample}` directory, but the main outputs of interest will be in `output/{sample}/mbased`. The outputs include:
1. The tabular results of the output `MBASED_expr_gene_results.txt`
2. The rds object of the MBASED raw output `MBASEDresults.rds`
3. The rds object of the MBASED gene output simplified dataframe `MBASED_expr_gene_results.txt`
4. Three descriptive figures of the results: `aseGenesBar.pdf`, `aseGenesDot.pdf`, `chromPlot.pdf` 
