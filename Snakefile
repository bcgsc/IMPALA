## Load config values
configfile: "config/defaults.yaml"
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"
configfile: "config/annotationPaths.yaml"

# path to the reference genome fasta
genome_path = config["genome"][config["genome_name"]]
genome_name = {"hg19": "GRCh37.75", "hg38": "GRCh38.100", "hg38_no_alt_TCGA_HTMCP_HPVs": "GRCh38.100"}[config["genome_name"]]

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
    log: "output/{sample}/log/phase_vcf_filter.log"
    shell:
        "zcat {input.phase} | grep -E '(PASS|#)' | grep -E '(0/1|\||#)' | awk '/^#/||length($4)==1 && length($5)==1' | bgzip > {output} 2> {log}"


rule phase_vcf_index:
    input:
        vcf = "output/{sample}/phase.het.pass.snps.vcf.gz"
    output:
        "output/{sample}/phase.het.pass.snps.vcf.gz.tbi"
    log: "output/{sample}/log/phase_vcf_index.log"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    shell:
        "tabix {input.vcf} &> {log}"

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
    log: "output/{sample}/log/rna_snv_calling.log"
    threads: 20
    shell:
        """
        configureStrelkaGermlineWorkflow.py --bam={input.bam} --referenceFasta={input.ref} --forcedGT={input.vcf} --rna --runDir=output/{wildcards.sample}/StrelkaRNA &> {log}
        output/{wildcards.sample}/StrelkaRNA/runWorkflow.py -m local -j {threads} &> {log}
        """

rule pass_filt:
    input:
        vcf="output/{sample}/StrelkaRNA/results/variants/genome.S1.vcf.gz"
    output:
        "output/{sample}/rna.forceGT.pass.vcf.gz"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    log: "output/{sample}/log/pass_filt.log"
    shell:
        "zcat {input.vcf} | grep -E '(PASS|#)' | bgzip > {output} 2> {log}"

rule rna_snv_index:
    input:
        vcf = "output/{sample}/rna.forceGT.pass.vcf.gz"
    output:
        "output/{sample}/rna.forceGT.pass.vcf.gz.tbi"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    log: "output/{sample}/log/rna_snv_index.log"
    shell:
        "tabix {input.vcf} &> {log}"

### -------------------------------------------------------------------
### Intersect RNA VCF with Phased VCF
### -------------------------------------------------------------------
rule intersect:
    input:
        vcf1 = "output/{sample}/phase.het.pass.snps.vcf.gz",
        vcf2 = "output/{sample}/rna.forceGT.pass.vcf.gz",
        index = "output/{sample}/rna.forceGT.pass.vcf.gz.tbi"
    output:
        "output/{sample}/rna.isec.snps.vcf"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/bcftools:1.15--h0ea216a_2"
    log: "output/{sample}/log/intersect.log"
    shell:
        """
        bcftools isec {input.vcf2} {input.vcf1} -p output/{wildcards.sample}/isec -n =2 -w 1 &> {log}
        mv output/{wildcards.sample}/isec/0000.vcf {output}
        """

rule intersect_gz:
    input:
        "output/{sample}/rna.isec.snps.vcf"
    output:
        "output/{sample}/rna.isec.snps.vcf.gz"
    conda: "config/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    log: "output/{sample}/log/intersect_gz.log"
    shell:
        "bgzip {input} &> {log}"

### -------------------------------------------------------------------
### Annotate and filter VCF with genes
### -------------------------------------------------------------------

rule snpEff:
    input: 
        "output/{sample}/rna.isec.snps.vcf"
    output:
        "output/{sample}/rna.isec.snps.snpEff.vcf"
    #singularity: "docker://quay.io/biocontainers/snpeff:5.1--hdfd78af_1"
    params:
        genome = genome_name,
        java = config["softwarePath"]["java"],
        snpEff = config["softwarePath"]["snpEff"],
        snpEff_config = config["annotationPath"]["snpEff_config"],
        snpEff_datadir = config["annotationPath"]["snpEff_datadir"]
    log: "output/{sample}/log/snpEff.log"
    shell:
        """
            {params.java} -Xmx64g -jar {params.snpEff} \
            -v {params.genome} \
            -c {params.snpEff_config} \
            -dataDir {params.snpEff_datadir} \
            -noStats \
            {input} > {output} 2> {log}
        """

rule dbSNP_annotation:
    input:
        "output/{sample}/rna.isec.snps.snpEff.vcf"
    output:
        "output/{sample}/rna.isec.snps.snpEff.dbSNP.vcf"
    params:
        java = config["softwarePath"]["java"],
        snpSift = config["softwarePath"]["snpSift"],
        dbSNP = config["annotationPath"]["dbSNP_database"]
    log: "output/{sample}/log/dbSNP_annotation.log"
    shell:
        """
        {params.java} -Xmx64g -jar {params.snpSift} annotate \
        {params.dbSNP} \
        {input} > {output} 2> {log}
        """

rule snpSift:
    input: 
        "output/{sample}/rna.isec.snps.snpEff.dbSNP.vcf"
    output:
        geneFilter = "output/{sample}/rna.isec.filterSnps.vcf",
        tsv = "output/{sample}/rna.isec.filterSnps.tsv"
    #singularity: "docker://quay.io/biocontainers/snpsift:4.2--hdfd78af_5"
    log: "output/{sample}/log/snpSift.log"
    params:
        java = config["softwarePath"]["java"],
        snpSift = config["softwarePath"]["snpSift"],
    shell:
        """
        {params.java} -Xmx64g -jar {params.snpSift} filter "( exists ANN[0].GENE )" {input} | \
        grep -E '(#|RS)' > {output.geneFilter} 2> {log}

        {params.java} -Xmx64g -jar {params.snpSift} extractFields {output.geneFilter} \
            CHROM POS GEN[0].AD ALT REF ANN[0].GENE ANN[0].BIOTYPE > {output.tsv} 2> {log}
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
        "bedtools intersect -loj -a {input.vcf} -b {input.bed} | cut -f 1-10,14,17,18 > {output} "

### -------------------------------------------------------------------
### Run MBASED
### -------------------------------------------------------------------

rule mbased:
    input:
        phase = lambda w: config["samples"][w.sample]["phase"],
        tsv = "output/{sample}/rna.isec.filterSnps.tsv",
    output:
        "output/{sample}/mBASED/MBASEDresults.rds"
    threads: 20
    log: "output/{sample}/log/mbased.log"
    shell:
        "scripts/mbased.snpEff.R --phase={input.phase} --rna={input.tsv} --outdir=output/{wildcards.sample}/mBASED &> {log}"

rule addExpression:
    input:
        rds = "output/{sample}/mBASED/MBASEDresults.rds",
        rpkm = rpkm_path
    output:
        "output/{sample}/mBASED/MBASED_expr_gene_results.txt"
    log: "output/{sample}/log/addExpression.log"
    shell:
        "scripts/addExpression.R --mbased={input.rds} --sample={wildcards.sample} --rpkm={input.rpkm} --min=1 --outdir=output/{wildcards.sample}/mBASED &> {log}"

rule figures:
    input:
        txt = "output/{sample}/mBASED/MBASED_expr_gene_results.txt",
        bed = gene_anno,
        rpkm = rpkm_path
    output:
        "output/{sample}/mBASED/sankeyPlot.html"
    log: "output/{sample}/log/figures.log"
    shell:
        "scripts/figures.R --mbased={input.txt} --rpkm={input.rpkm} --gene={input.bed} --sample={wildcards.sample} --outdir=output/{wildcards.sample}/mBASED &> {log}"
