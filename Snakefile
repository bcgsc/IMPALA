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

# Check phased
phased = config["phased"]
### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
    input:
	    expand("output/{sample}/mBASED/chromPlot.pdf",sample=sample_ids)


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

if phased:
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
        	configureStrelkaGermlineWorkflow.py \
			--bam={input.bam} \
			--referenceFasta={input.ref} \
			--forcedGT={input.vcf} \
			--rna \
			--runDir=output/{wildcards.sample}/StrelkaRNA &> {log}
        	output/{wildcards.sample}/StrelkaRNA/runWorkflow.py -m local -j {threads} &> {log}
        	"""
else:
	rule rna_snv_calling:
            input:
                bam = lambda w: config["samples"][w.sample]["rna"],
                ref = genome_path,
        output:
                "output/{sample}/StrelkaRNA/results/variants/genome.S1.vcf.gz"
        conda: "config/strelka.yaml"
        singularity: "docker://quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"
        log: "output/{sample}/log/rna_snv_calling.log"
        threads: 20
        shell:
                """
                configureStrelkaGermlineWorkflow.py \
                        --bam={input.bam} \
                        --referenceFasta={input.ref} \
                        --rna \
                        --runDir=output/{wildcards.sample}/StrelkaRNA &> {log}
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
if phased:
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
		input: "output/{sample}/rna.isec.snps.vcf"
    		output: "output/{sample}/rna.isec.snps.vcf.gz"
    		conda: "config/ase-env.yaml"
    		singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    		log: "output/{sample}/log/intersect_gz.log"
    		shell: "bgzip {input} &> {log}"
	vcf = "output/{sample}/rna.isec.snps.vcf.gz"
else:
	vcf = "output/{sample}/rna.forceGT.pass.vcf.gz"


### -------------------------------------------------------------------
### Annotate and filter VCF with genes
### -------------------------------------------------------------------

rule snpEff:
    input: vcf
    output:
        "output/{sample}/rna.isec.snps.snpEff.vcf"
    singularity: "docker://quay.io/biocontainers/snpeff:5.0--hdfd78af_1"
    params:
        genome = genome_name,
        snpEff_config = config["annotationPath"]["snpEff_config"],
        snpEff_datadir = config["annotationPath"]["snpEff_datadir"]
    log: "output/{sample}/log/snpEff.log"
    shell:
        """
            snpEff -Xmx64g \
            -v {params.genome} \
            -c {params.snpEff_config} \
            -dataDir {params.snpEff_datadir} \
            -noStats \
            {input} > {output} 2> {log}
        """

rule snpSift:
    input: "output/{sample}/rna.isec.snps.snpEff.vcf"
    output:
        geneFilter = "output/{sample}/rna.isec.filterSnps.vcf",
        tsv = "output/{sample}/rna.isec.filterSnps.tsv"
    singularity: "docker://quay.io/biocontainers/snpsift:5.1d--hdfd78af_0"
    log: "output/{sample}/log/snpSift.log"
    shell:
        """
        SnpSift -Xmx64g filter "( exists ANN[0].GENE )" {input} > {output.geneFilter} 2> {log}

        SnpSift -Xmx64g extractFields {output.geneFilter} \
            CHROM POS GEN[0].AD ALT REF ANN[0].GENE ANN[0].BIOTYPE > {output.tsv} 2> {log}
        """

### -------------------------------------------------------------------
### Run MBASED
### -------------------------------------------------------------------

if phased:
	rule mbased:
    		input:
        		phase = lambda w: config["samples"][w.sample]["phase"],
        		tsv = "output/{sample}/rna.isec.filterSnps.tsv",
    		output:
        		"output/{sample}/mBASED/MBASEDresults.rds"
    		threads: 20
    		log: "output/{sample}/log/mbased.log"
    		shell:
        		"""
			scripts/mbased.snpEff.R \
				--threads={threads} \
				--phase={input.phase} \
				--rna={input.tsv} \
				--outdir=output/{wildcards.sample}/mBASED &> {log}
			"""
else:
	rule mbased:
                input:
                        tsv = "output/{sample}/rna.isec.filterSnps.tsv",
                output:
                        "output/{sample}/mBASED/MBASEDresults.rds"
                threads: 20
                log: "output/{sample}/log/mbased.log"
                shell:
                        """
                        scripts/mbased.snpEff.R \
                                --threads={threads} \
                                --rna={input.tsv} \
                                --outdir=output/{wildcards.sample}/mBASED &> {log}
                        """
	

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
        "output/{sample}/mBASED/chromPlot.pdf"
    log: "output/{sample}/log/figures.log"
    shell:
        "scripts/figures.R --mbased={input.txt} --rpkm={input.rpkm} --gene={input.bed} --sample={wildcards.sample} --outdir=output/{wildcards.sample}/mBASED &> {log}"
