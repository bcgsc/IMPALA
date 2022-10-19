## Load config values
configfile: "config/defaults.yaml"
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"

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

# Check cancer analysis
cancer_analysis = config["cancer_analysis"]

if cancer_analysis:
	ruleorder: summaryTableCancer > summaryTable
else:
	ruleorder: summaryTable > summaryTableCancer

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
		expand("output/{sample}/figures/chromPlot.pdf",sample=sample_ids),
		expand("output/{sample}/summaryTable.tsv", sample=sample_ids)


### -------------------------------------------------------------------
### Alignment and expression matrix
### -------------------------------------------------------------------

if phased:
	rule starAlignment:
		input:
			vcf = lambda w: config["samples"][w.sample]["phase"],
			r1 = lambda w: config["samples"][w.sample]["R1"],
			r2 = lambda w: config["samples"][w.sample]["R2"]
		output: 
			"output/{sample}/0_alignment/starAligned.sortedByCoord.out.bam",
			"output/{sample}/0_alignment/starAligned.toTranscriptome.out.bam"
		singularity: "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
		threads: config["threads"]
		params: 
			ref = config["starReferencePath"]
		log: "output/{sample}/log/starAlignment.log"
		shell:
			"""
			mkdir -p output/{wildcards.sample}/0_alignment
			STAR \
				--genomeDir {params.ref} \
				--runThreadN {threads} \
				--readFilesIn  {input.r1} {input.r2} \
				--outFileNamePrefix output/{wildcards.sample}/0_alignment/star \
				--outSAMtype BAM SortedByCoordinate \
				--outSAMunmapped Within \
				--outSAMattributes Standard \
				--waspOutputMode SAMtag \
				--varVCFfile {input.vcf} \
				--quantMode TranscriptomeSAM \
				--twopassMode Basic \
				--twopass1readsN -1 &> {log}
			"""
	
	rule waspFilter:
		input: "output/{sample}/0_alignment/starAligned.sortedByCoord.out.bam"
		output: "output/{sample}/0_alignment/starAligned.waspFilter.bam"
		singularity: "docker://quay.io/biocontainers/samtools:1.16.1--h6899075_1"
		log: "output/{sample}/log/waspFilter.log"
		shell:
			"""
			samtools view -h {input} | grep -e '^@' -e 'vW:i:1' | samtools view -b -S > {output} 2> {log}
			samtools index {output} &> {log}
			"""
	alignmentFile = "output/{sample}/0_alignment/starAligned.waspFilter.bam"

else:
	rule starAlignment:
		input:
			r1 = lambda w: config["samples"][w.sample]["R1"],
			r2 = lambda w: config["samples"][w.sample]["R2"]
		output: 
			"output/{sample}/0_alignment/starAligned.sortedByCoord.out.bam",
			"output/{sample}/0_alignment/starAligned.toTranscriptome.out.bam"
		singularity: "docker://quay.io/biocontainers/star:2.7.10a--h9ee0642_0"
		threads: config["threads"]
		params: 
			ref = config["starReferencePath"]
		log: "output/{sample}/log/starAlignment.log"
		shell:
			"""
			mkdir -p output/{wildcards.sample}/0_alignment
			STAR \
				--genomeDir {params.ref} \
				--runThreadN {threads} \
				--readFilesIn  {input.r1} {input.r2} \
				--outFileNamePrefix output/{wildcards.sample}/0_alignment/star \
				--outSAMtype BAM SortedByCoordinate \
				--outSAMunmapped Within \
				--outSAMattributes Standard \
				--quantMode TranscriptomeSAM \
				--twopassMode Basic \
				--twopass1readsN -1 &> {log}
			"""
	alignmentFile = "output/{sample}/0_alignment/star/starAligned.sortedByCoord.out.bam"

rule rsem:
	input: "output/{sample}/0_alignment/starAligned.toTranscriptome.out.bam"
	output: "output/{sample}/0_alignment/star.genes.results"
	params:
			ref=config["rsemReferencePath"]
	singularity: "docker://quay.io/biocontainers/rsem:1.3.3--pl5321hecb563c_4"
	threads: config["threads"]
	log: "output/{sample}/log/rsem.log"
	shell:
			"""
			rsem-calculate-expression \
					--alignments \
					-p {threads} \
					--paired-end \
					{input} \
					{params.ref} \
					output/{wildcards.sample}/0_alignment/star &> {log}
			"""
			
rule rsemExpressionMatrix:
	input: "output/{sample}/0_alignment/star.genes.results"
	output: "output/{sample}/0_alignment/expression_matrix.tsv"
	params:
		annotation = "annotation/biomart_ensembl100_GRCh38.sorted.bed"
	singularity: "docker://glenn032787/ase_rcontainer:1.0"
	log: "output/{sample}/log/rsemExpressionMatrix.log"
	shell:
		"""
		Rscript scripts/generateExpressionMatrix.R \
			-i {input} \
			-o output/{wildcards.sample}/0_alignment \
			-a {params.annotation} \
			-s {wildcards.sample} &> {log}
		"""
### -------------------------------------------------------------------
### Call and filter the phase vcf
### -------------------------------------------------------------------

rule phase_vcf_filter:
    input: 
        phase = lambda w: config["samples"][w.sample]["phase"]
    output:
        "output/{sample}/1_variant/phase.het.pass.snps.vcf.gz"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    log: "output/{sample}/log/phase_vcf_filter.log"
    shell:
        """
		mkdir -p output/{wildcards.sample}/1_variant
		cat {input.phase} | grep -E '(PASS|#)' | grep -E '(0/1|\||#)' | awk '/^#/||length($4)==1 && length($5)==1' | bgzip > {output} 2> {log}
		"""


rule phase_vcf_index:
    input:
        vcf = "output/{sample}/1_variant/phase.het.pass.snps.vcf.gz"
    output:
        "output/{sample}/1_variant/phase.het.pass.snps.vcf.gz.tbi"
    log: "output/{sample}/log/phase_vcf_index.log"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    shell:
        "tabix {input.vcf} &> {log}"

### -------------------------------------------------------------------
### Call and filter the RNA SNVs
### -------------------------------------------------------------------

def getAlignment(wildcards):
	if config["samples"][wildcards.sample]["rna"] == None:
		print(alignmentFile)
		return alignmentFile
	return config["samples"][wildcards.sample]["rna"]


if phased:
	rule rna_snv_calling:
		input:
			#bam = lambda w: config["samples"][w.sample]["rna"],
			bam = getAlignment,
        	vcf = "output/{sample}/1_variant/phase.het.pass.snps.vcf.gz",
        	ref = genome_path,
			index = "output/{sample}/1_variant/phase.het.pass.snps.vcf.gz.tbi"
		output:
			temp("output/{sample}/StrelkaRNA/results/variants/genome.S1.vcf.gz")
		conda: "config/conda/strelka.yaml"
		singularity: "docker://quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"
		log: "output/{sample}/log/rna_snv_calling.log"
		threads: config["threads"]
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
                temp("output/{sample}/StrelkaRNA/results/variants/genome.S1.vcf.gz")
        conda: "config/conda/strelka.yaml"
        singularity: "docker://quay.io/biocontainers/strelka:2.9.10--h9ee0642_1"
        log: "output/{sample}/log/rna_snv_calling.log"
        threads: config["threads"]
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
        temp("output/{sample}/1_variant/rna.forceGT.pass.vcf.gz")
    conda: "config/conda/ase-env.yaml"
    singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    log: "output/{sample}/log/pass_filt.log"
    shell:
        "zcat {input.vcf} | grep -E '(PASS|#)' | bgzip > {output} 2> {log} && rm -rf output/{wildcards.sample}/StrelkaRNA/"

rule rna_snv_index:
    input:
        vcf = "output/{sample}/1_variant/rna.forceGT.pass.vcf.gz"
    output:
        temp("output/{sample}/1_variant/rna.forceGT.pass.vcf.gz.tbi")
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
        		vcf1 = "output/{sample}/1_variant/phase.het.pass.snps.vcf.gz",
        		vcf2 = "output/{sample}/1_variant/rna.forceGT.pass.vcf.gz",
        		index = "output/{sample}/1_variant/rna.forceGT.pass.vcf.gz.tbi"
    		output:
        		temp("output/{sample}/1_variant/rna.isec.snps.vcf")
    		conda: "config/conda/ase-env.yaml"
    		singularity: "docker://quay.io/biocontainers/bcftools:1.15--h0ea216a_2"
    		log: "output/{sample}/log/intersect.log"
    		shell:
        		"""
        		bcftools isec {input.vcf2} {input.vcf1} -p output/{wildcards.sample}/isec -n =2 -w 1 &> {log} 
        		mv output/{wildcards.sample}/isec/0000.vcf {output} 
				rm -rf output/{wildcards.sample}/isec
        		"""

	rule intersect_gz:
		input: "output/{sample}/1_variant/rna.isec.snps.vcf"
    		output: temp("output/{sample}/1_variant/rna.isec.snps.vcf.gz")
    		conda: "config/conda/ase-env.yaml"
    		singularity: "docker://quay.io/biocontainers/htslib:1.15--h9753748_0"
    		log: "output/{sample}/log/intersect_gz.log"
    		shell: "bgzip {input} &> {log}"
	vcf = "output/{sample}/1_variant/rna.isec.snps.vcf.gz"
else:
	vcf = "output/{sample}/1_variant/rna.forceGT.pass.vcf.gz"


### -------------------------------------------------------------------
### Annotate and filter VCF with genes
### -------------------------------------------------------------------

rule snpEff:
    input: vcf
    output:
        temp("output/{sample}/1_variant/rna.isec.snps.snpEff.vcf")
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
    input: "output/{sample}/1_variant/rna.isec.snps.snpEff.vcf"
    output:
        geneFilter = "output/{sample}/1_variant/rna.isec.filterSnps.vcf",
        tsv = "output/{sample}/1_variant/rna.isec.filterSnps.tsv"
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
        		tsv = "output/{sample}/1_variant/rna.isec.filterSnps.tsv"
    		output:
        		"output/{sample}/2_mBASED/MBASEDresults.rds"
    		threads: config["threads"]
    		log: "output/{sample}/log/mbased.log"
		singularity: "docker://glenn032787/ase_rcontainer:1.0"
    		shell:
        		"""
				Rscript scripts/mbased.snpEff.R \
					--threads={threads} \
					--phase={input.phase} \
					--rna={input.tsv} \
					--outdir=output/{wildcards.sample}/2_mBASED &> {log}
			"""
else:
	rule mbased:
			input:
				tsv = "output/{sample}/1_variant/rna.isec.filterSnps.tsv",
			output:
				"output/{sample}/2_mBASED/MBASEDresults.rds"
			threads: config["threads"]
			singularity: "docker://glenn032787/ase_rcontainer:1.0"
			log: "output/{sample}/log/mbased.log"
			shell:
				"""
				Rscript scripts/mbased.snpEff.R \
						--threads={threads} \
						--rna={input.tsv} \
						--outdir=output/{wildcards.sample}/2_mBASED &> {log}
				"""

def getExpressionMatrix(wildcards):
	if config["samples"][wildcards.sample]["rna"] == None:
		return "output/{sample}/0_alignment/expression_matrix.tsv"
	return rpkm_path

rule addExpression:
	input:
		rds = "output/{sample}/2_mBASED/MBASEDresults.rds",
		rpkm = getExpressionMatrix
	output: "output/{sample}/2_mBASED/MBASED_expr_gene_results.txt"
	singularity: "docker://glenn032787/ase_rcontainer:1.0"
	params:
		maf = config["maf_threshold"]
	log: "output/{sample}/log/addExpression.log"
	shell:
		"""
		Rscript scripts/addExpression.R \
			--mbased={input.rds} \
			--sample={wildcards.sample} \
			--rpkm={input.rpkm} \
			--min=1 \
			--maf_threshold={params.maf} \
			--outdir=output/{wildcards.sample}/2_mBASED &> {log}
		"""

rule figures:
	input:
		txt = "output/{sample}/2_mBASED/MBASED_expr_gene_results.txt",
		bed = gene_anno,
		rpkm = getExpressionMatrix
	output:
		"output/{sample}/figures/chromPlot.pdf"
	singularity: "docker://glenn032787/ase_rcontainer:1.0"
	params:
		maf = config["maf_threshold"]
	log: "output/{sample}/log/figures.log"
	shell:
		"""
		mkdir -p output/{wildcards.sample}/figures

		Rscript scripts/figures.R \
			--mbased={input.txt} \
			--rpkm={input.rpkm} \
			--gene={input.bed} \
			--sample={wildcards.sample} \
			--maf_threshold={params.maf} \
			--outdir=output/{wildcards.sample}/figures &> {log}
		"""

rule summaryTable:
	input: "output/{sample}/2_mBASED/MBASED_expr_gene_results.txt"
	output: "output/{sample}/summaryTable.tsv"
	shell:
		"""
		cut -f 1,3,4,5,8,11 {input} > {output}
		"""


### -------------------------------------------------------------------
### Cancer analysis
### -------------------------------------------------------------------

rule annotateGenes:
	input: "output/{sample}/2_mBASED/MBASED_expr_gene_results.txt"
	output: 
		bed = "output/{sample}/3_cancer/raw/gene_annotation.bed",
		gene = "output/{sample}/3_cancer/raw/phasedGenes.txt"
	params:
		annotation = "annotation/biomart_ensembl100_GRCh38.sorted.bed"
	log: "output/{sample}/log/annotateGenes.log"
	shell:	
		"""
		mkdir -p output/{wildcards.sample}/3_cancer/raw

		cat {input} | cut -f1 > {output.gene} 2> {log}
		awk 'NR == FNR {{ keywords[$1]=1; next; }} {{ if ($4 in keywords) print; }}' {output.gene} {params.annotation} > {output.bed} 2> {log}
		"""


rule genomeLength:
	input: genome_path + ".fai"
	output: temp("output/{sample}/3_cancer/raw/genome.length")
	log: "output/{sample}/log/genomeLength.log"
	shell:
		"""
		cut -f 1,2 {input} > {output} 2> {log}
		"""
	


rule promoterSlop:
	input: 
		gene = "output/{sample}/3_cancer/raw/gene_annotation.bed",
		length = "output/{sample}/3_cancer/raw/genome.length"
	output: "output/{sample}/3_cancer/raw/promoter_annotation.bed"
	singularity: "docker://quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
	log: "output/{sample}/log/promoterSlop.log"
	shell:
		"""
		bedtools slop -l 2000 -r 500 -i {input.gene} -g {input.length} > {output} 2> {log}
		"""

rule cnv_dmr:
	input:
		cnv = lambda w: config["samples"][w.sample]["cnv"],
		methyl = lambda w: config["samples"][w.sample]["methyl"] 
	output: 
		cnv = "output/{sample}/3_cancer/raw/cnv.bed",
		methyl = "output/{sample}/3_cancer/raw/methyl.bed"
	singularity: "docker://glenn032787/ase_rcontainer:1.0"
	log: "output/{sample}/log/cnv_dmr.log"
	shell:
		"""
		mkdir -p output/{wildcards.sample}/3_cancer/raw 

		Rscript scripts/cnv_dmr_process.R \
			--methyl={input.methyl} \
			--cnv={input.cnv} \
			--outdir=output/{wildcards.sample}/3_cancer/raw &> {log}
		"""	

rule cnvIntersect:
	input: 
		cnv = "output/{sample}/3_cancer/raw/cnv.bed",
		gene = "output/{sample}/3_cancer/raw/gene_annotation.bed"
	output: "output/{sample}/3_cancer/intersect/cnv_intersect.bed"
	singularity: "docker://quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
	log: "output/{sample}/log/cnvIntersect.log"
	shell:
		"""
		mkdir -p output/{wildcards.sample}/3_cancer/intersect
		bedtools intersect -loj -a {input.gene} -b {input.cnv} | awk '$9 != "." {{print $0}}' > {output} 2> {log}
		"""


rule methylIntersect:
	input: 
		methyl = "output/{sample}/3_cancer/raw/methyl.bed",
		gene = "output/{sample}/3_cancer/raw/gene_annotation.bed"
	output: "output/{sample}/3_cancer/intersect/methyl_intersect.bed"
	singularity: "docker://quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
	log: "output/{sample}/log/methylIntersect.log"
	shell:
		"""
		mkdir -p output/{wildcards.sample}/3_cancer/intersect
		bedtools intersect -loj -a {input.gene} -b {input.methyl} | awk '$9 != "." {{print $0}}' > {output} 2> {log}
		"""

rule summaryTableCancer:
	input:
		cnv = "output/{sample}/3_cancer/intersect/cnv_intersect.bed",
		methyl = "output/{sample}/3_cancer/intersect/methyl_intersect.bed", 
		ase = "output/{sample}/2_mBASED/MBASED_expr_gene_results.txt" 
	output: "output/{sample}/summaryTable.tsv"
	singularity: "docker://glenn032787/ase_rcontainer:1.0"
	log: "output/{sample}/log/summaryTable.log"
	shell:
		"""
		Rscript scripts/summaryTable.R \
			--cnv={input.cnv} \
			--methyl={input.methyl} \
			--ase={input.ase} \
			--sample={wildcards.sample} \
			--outdir=output/{wildcards.sample} 2> {log}
		"""


