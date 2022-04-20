## Load config values
configfile: "config/defaults.yaml"
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"

# path to the reference genome fasta
genome_path = config["genome"][config["genome_name"]]
genome_name = {"hg19": "GRCh37.75", "hg38": "GRCh38.99", "hg38_no_alt_TCGA_HTMCP_HPVs": "GRCh38.99"}[config["genome_name"]]

# Ensembl 100 genes
gene_anno = config["annotation"][config["genome_name"]]

# path to the RPKM matrix
rpkm_path = config["matrix"]

# samples
samples_dict = config["samples"]
sample_ids = samples_dict.keys()

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
    input:
	    expand("output/{sample}/mBASED/sankeyPlot.html",sample=sample_ids)

### -------------------------------------------------------------------
### Call and filter the DNA SNVs # Unused rules (Used phased VCF instead)
### -------------------------------------------------------------------

rule dna_snv_calling: # UNUSED RULE
    input:
        bam = lambda w: config["samples"][w.sample]["dna"],
        ref = genome_path,
        bed = gene_anno
    output:
        "output/{sample}/StrelkaDNA/results/variants/genome.S1.vcf.gz"
    conda: "config/strelka.yaml"
    singularity: "docker://quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"
    threads: 20
    shell:
        """
        configureStrelkaGermlineWorkflow.py --bam={input.bam} --referenceFasta={input.ref} --callRegions={input.bed} --rna --runDir=output/{wildcards.sample}/StrelkaDNA
        output/{wildcards.sample}/StrelkaDNA/runWorkflow.py -m local -j {threads}
        """

rule dna_snv_filt: # UNUSED RULE
    input:
        vcf = "output/{sample}/StrelkaDNA/results/variants/genome.S1.vcf.gz"
    output:
        "output/{sample}/dna.het.pass.snps.vcf.gz"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    shell:
        "zcat {input.vcf} | grep -E '(PASS|#)' | grep -E '(0/1|#)' | awk '/^#/||length($4)==1 && length($5)==1' | bgzip > {output}"

rule dna_snv_index: # UNUSED RULE
    input:
        vcf = "output/{sample}/dna.het.pass.snps.vcf.gz"
    output:
        "output/{sample}/dna.het.pass.snps.vcf.gz.tbi"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    shell:
        "tabix {input.vcf}"

### -------------------------------------------------------------------
### Call and filter the phase vcf
### -------------------------------------------------------------------
rule phase_vcf_filter:
    input: 
        phase = lambda w: config["samples"][w.sample]["phase"]
    output:
        "output/{sample}/phase.het.pass.snps.vcf.gz"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    shell:
        "cat {input.phase} | grep -E '(PASS|#)' | grep -E '(0/1|\||#)' | awk '/^#/||length($4)==1 && length($5)==1' | bgzip > {output}"


rule phase_vcf_index:
    input:
        vcf = "output/{sample}/phase.het.pass.snps.vcf.gz"
    output:
        "output/{sample}/phase.het.pass.snps.vcf.gz.tbi"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    shell:
        "tabix {input.vcf}"

### -------------------------------------------------------------------
### Call and filter the RNA SNVs
### -------------------------------------------------------------------

rule rna_snv_calling:
    input:
        bam = lambda w: config["samples"][w.sample]["rna"],
        vcf = "output/{sample}/phase.het.pass.snps.vcf.gz",
        ref = genome_path,
        index = "output/{sample}/phase.het.pass.snps.vcf.gz.tbi"
    output:
        "output/{sample}/StrelkaRNA/results/variants/genome.S1.vcf.gz"
    conda: "config/strelka.yaml"
    singularity: "docker://quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"
    threads: 20
    shell:
        """
        configureStrelkaGermlineWorkflow.py --bam={input.bam} --referenceFasta={input.ref} --forcedGT={input.vcf} --rna --runDir=output/{wildcards.sample}/StrelkaRNA
        output/{wildcards.sample}/StrelkaRNA/runWorkflow.py -m local -j {threads}
        """

rule pass_filt:
    input:
        vcf="output/{sample}/StrelkaRNA/results/variants/genome.S1.vcf.gz"
    output:
        "output/{sample}/rna.forceGT.pass.vcf.gz"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    shell:
        "zcat {input.vcf} | grep -E '(PASS|#)' | bgzip > {output}"

rule rna_snv_index:
    input:
        vcf = "output/{sample}/rna.forceGT.pass.vcf.gz"
    output:
        "output/{sample}/rna.forceGT.pass.vcf.gz.tbi"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    shell:
        "tabix {input.vcf}"

### -------------------------------------------------------------------
### Intersect RNA VCF with Phased VCF
### -------------------------------------------------------------------
rule intersect:
    input:
        vcf1 = "output/{sample}/phase.het.pass.snps.vcf.gz",
        vcf2 = "output/{sample}/rna.forceGT.pass.vcf.gz",
        index = "output/{sample}/rna.forceGT.pass.vcf.gz.tbi"
    output:
        "output/{sample}/rna.isec.dna.snps.vcf"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/bcftools:1.15--h0ea216a_2"
    shell:
        """
        bcftools isec {input.vcf2} {input.vcf1} -p output/{wildcards.sample}/isec -n =2 -w 1
        mv output/{wildcards.sample}/isec/0000.vcf {output}
        """

rule intersect_gz:
    input:
        "output/{sample}/rna.isec.dna.snps.vcf"
    output:
        "output/{sample}/rna.isec.dna.snps.vcf.gz"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    shell:
        "bgzip {input}"

### -------------------------------------------------------------------
### Annotate and filter VCF with genes
### -------------------------------------------------------------------

rule snpEff:
    input: 
        "output/{sample}/rna.isec.dna.snps.vcf"
    output:
        "output/{sample}/rna.isec.dna.snps.annot.vcf"
    singularity: "docker://quay.io/biocontainers/snpeff:5.1--hdfd78af_1"
    params:
        genome = genome_name
    shell:
        """
        java -Xmx100g -jar /usr/local/share/snpeff-5.1-1/snpEff.jar {params.genome} {input} > {output}
        """

rule snpSift:
    input: 
        "output/{sample}/rna.isec.dna.snps.annot.vcf"
    output:
        geneFilter = "output/{sample}/rna.isec.dna.snps.annotGene.vcf",
        tsv = "output/{sample}/rna.isec.dna.snps.annotGene.tsv"
    singularity: "docker://quay.io/biocontainers/snpsift:4.2--hdfd78af_5"
    shell:
        """
        java -Xmx100g -jar /usr/local/share/snpsift-4.3.1t-3/SnpSift.jar filter "( ANN[0].GENE exist )" {input} > {output.geneFilter}
        java -Xmx100g -jar /usr/local/share/snpsift-4.3.1t-3/SnpSift.jar extractFields {output.geneFilter} CHROM POS GEN[0].AD ALT REF ANN[0].GENE ANN[0].BIOTYPE > {output.tsv}
        """

rule intersect_genes: # UNUSED RULE
    input:
        vcf = "output/{sample}/rna.isec.dna.snps.vcf.gz",
        bed = gene_anno
    output:
        "output/{sample}/rna.isec.dna.snps.genes.vcf.gz"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
    shell:
        "bedtools intersect -loj -a {input.vcf} -b {input.bed} | cut -f 1-10,14,17,18 > {output}"

### -------------------------------------------------------------------
### Run MBASED
### -------------------------------------------------------------------

rule mbased:
    input:
        phase = lambda w: config["samples"][w.sample]["phase"],
        tsv = "output/{sample}/rna.isec.dna.snps.genes.vcf.gz",
    output:
        "output/{sample}/mBASED/MBASEDresults.rds"
    threads: 20
    shell:
        "scripts/mbased.R --phase={input.phase} --rna={input.tsv} --outdir=output/{wildcards.sample}/mBASED"

rule addExpression:
    input:
        rds = "output/{sample}/mBASED/MBASEDresults.rds",
        rpkm = rpkm_path
    output:
        "output/{sample}/mBASED/MBASED_expr_gene_results.txt"
    shell:
        "scripts/addExpression.R --mbased={input.rds} --sample={wildcards.sample} --rpkm={input.rpkm} --min=1 --outdir=output/{wildcards.sample}/mBASED"

rule figures:
    input:
        txt = "output/{sample}/mBASED/MBASED_expr_gene_results.txt",
        bed = gene_anno,
        rpkm = rpkm_path
    output:
        "output/{sample}/mBASED/sankeyPlot.html"
    shell:
        "scripts/figures.R --mbased={input.txt} --rpkm={input.rpkm} --gene={input.bed} --sample={wildcards.sample} --outdir=output/{wildcards.sample}/mBASED"
