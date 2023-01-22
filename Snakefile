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

# Possible tissues
possibleTissue = ["Adipose", "Adrenal_gland", "allTissue", "Bladder", "Blood", "Brain", "Breast", 
		  "Cervix", "Colon", "Esophagus", "Fallopian_tube", "Heart", "Kidney", "Liver",
		  "Lung", "Muscle", "Nerve", "Ovary", "Pancreas", "Pituitary", "Prostate", "Salivary",
		  "Skin", "Small_intestine", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus","Vagina"]

# Check phased
phased = config["phased"]

# Check cancer analysis
cancer_analysis = config["cancer_analysis"]

if cancer_analysis:
	figure = expand("output/{sample}/figures/aseCause.pdf", sample=sample_ids)
else:
	figure = []

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
		expand("output/{sample}/params.txt",sample=sample_ids),
		expand("output/{sample}/figures/sankeyPlot.html",sample=sample_ids),
		expand("output/{sample}/summaryTable.tsv", sample=sample_ids),
		expand("output/{sample}/figures/karyogram.pdf", sample = sample_ids),
		figure

### -------------------------------------------------------------------
### Params
### -------------------------------------------------------------------

rule params:
	output: "output/{sample}/params.txt"
	params: 
		sample_info = lambda w: config["samples"][w.sample],
		genome_name = config["genome_name"],
		expressionMatrix = config["matrix"],
		phase = config["phased"],
		cancer = config["cancer_analysis"],
		maf = config["maf_threshold"],
		threads = config["threads"]
	shell:
		"""
		date >> {output}
		printf "{wildcards.sample}\n------------------------\n\n" >> {output}
		echo {params.sample_info}| sed 's/{{//' | sed 's/}}//' | tr  , '\n'  >> {output}
		echo expression_matrix: {params.expressionMatrix} >> {output}
		printf "\n------------------------\n" >> {output}
		echo Major Allele Frequency Threshold: {params.maf} >> {output}
		echo Genome: {params.genome_name} >> {output}
		echo Phased: {params.phase}  >> {output}
		echo Cancer analysis: {params.cancer} >> {output}	
		echo Threads: {params.threads} >> {output}
		"""

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
	singularity: "docker://glenn032787/ase_rcontainer:2.0"
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
		zcat {input.phase} | grep -E '(PASS|#)' | grep -E '(0/1|\||#)' | awk '/^#/||length($4)==1 && length($5)==1' | bgzip > {output} 2> {log}
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
	if "rna" not in config["samples"][wildcards.sample] or config["samples"][wildcards.sample]["rna"] == None:
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
    #conda: "config/conda/ase-env.yaml"
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
    		#conda: "config/conda/ase-env.yaml"
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
    		#conda: "config/conda/ase-env.yaml"
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
        genome = config["annotationPath"]["snpEff_genomeName"],
        snpEff_config = config["annotationPath"]["snpEff_config"],
        snpEff_datadir = config["annotationPath"]["snpEff_datadir"]
    log: "output/{sample}/log/snpEff.log"
    shell:
        """
        snpEff -Xmx16g \
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
        SnpSift -Xmx16g filter "( exists ANN[0].GENE )" {input} > {output.geneFilter} 2> {log}

        SnpSift -Xmx16g extractFields {output.geneFilter} \
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
		singularity: "docker://glenn032787/ase_rcontainer:2.0"
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
			singularity: "docker://glenn032787/ase_rcontainer:2.0"
			log: "output/{sample}/log/mbased.log"
			shell:
				"""
				Rscript scripts/mbased.snpEff.R \
						--threads={threads} \
						--rna={input.tsv} \
						--outdir=output/{wildcards.sample}/2_mBASED &> {log}
				"""

def getExpressionMatrix(wildcards):
	if "rna" not in config["samples"][wildcards.sample] or  config["samples"][wildcards.sample]["rna"] == None:
		return "output/{sample}/0_alignment/expression_matrix.tsv"
	return rpkm_path

rule addExpression:
	input:
		rds = "output/{sample}/2_mBASED/MBASEDresults.rds",
		rpkm = getExpressionMatrix
	output: "output/{sample}/2_mBASED/MBASED_expr_gene_results.txt"
	singularity: "docker://glenn032787/ase_rcontainer:2.0"
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
		"output/{sample}/figures/sankeyPlot.html"
	singularity: "docker://glenn032787/ase_rcontainer:2.0"
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



### -------------------------------------------------------------------
### Cancer analysis - Copy Number Variant
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
		awk 'NR == FNR {{ keywords[$1]=1; next; }} {{ if ($4 in keywords) print; }}' {output.gene} {params.annotation} | \
			cut -f 1,2,3,4 | uniq > {output.bed} 2> {log}
		"""


rule genomeLength:
	input: genome_path + ".fai"
	output: temp("output/{sample}/3_cancer/raw/genome.length")
	log: "output/{sample}/log/genomeLength.log"
	shell:
		"""
		cut -f 1,2 {input} > {output} 2> {log}
		"""
	


rule promoterFlank:
	input: 
		gene = "output/{sample}/3_cancer/raw/gene_annotation.bed",
		length = "output/{sample}/3_cancer/raw/genome.length"
	output: "output/{sample}/3_cancer/raw/promoter_annotation.bed"
	singularity: "docker://quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
	log: "output/{sample}/log/promoterFlank.log"
	shell:
		"""
		bedtools flank -l 2000 -r 500 -i {input.gene} -g {input.length} > {output} 2> {log}
		"""

def getTumorContent(wildcards):
	if "tumorContent" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["tumorContent"] != None:
 		return config["samples"][wildcards.sample]["tumorContent"]
	else:
		return 1.0

rule cnv_preprocess:
	input:  lambda w: config["samples"][w.sample]["cnv"]
	output: "output/{sample}/3_cancer/raw/cnv.bed"
	params: 
		tumor = getTumorContent
	singularity: "docker://glenn032787/ase_rcontainer:2.0"
	log: "output/{sample}/log/cnv_preprocess.log"
	shell:
		"""
		Rscript scripts/cnv_preprocess.R \
			--cnv={input} \
			--tumorContent={params.tumor} \
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
		bedtools intersect -loj -a {input.gene} -b {input.cnv} | awk '$10 != "." {{print $0}}' > {output} 2> {log}
		"""


### -------------------------------------------------------------------
### Cancer analysis - Methylation
### -------------------------------------------------------------------

rule dmr_preprocess:
        input: lambda w: config["samples"][w.sample]["methyl"]
        output: "output/{sample}/3_cancer/raw/methyl.bed"
        singularity: "docker://glenn032787/ase_rcontainer:2.0"
        log: "output/{sample}/log/dmr_preprocess.log"
        shell:
                """
                Rscript scripts/methyl_preprocess.R \
                        --methyl={input} \
                        --outdir=output/{wildcards.sample}/3_cancer/raw &> {log}
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
		bedtools intersect -loj -a {input.gene} -b {input.methyl} | awk '$10 != "." {{print $0}}' > {output} 2> {log}
		"""



### -------------------------------------------------------------------
### Cancer analysis - Transciption Factor Binding Site
### -------------------------------------------------------------------

rule promoterID:
	input: "output/{sample}/3_cancer/raw/promoter_annotation.bed"
	output:
		id2gene = temp("output/{sample}/3_cancer/tfbs/id2gene.txt"),
		id = temp("output/{sample}/3_cancer/tfbs/id.txt")
	log: "output/{sample}/log/promoterID.log"
	shell:
		"""
		mkdir -p output/{wildcards.sample}/3_cancer/tfbs
		awk 'BEGIN {{print "sequence_id\tgene"}}; {{print $1 ":" $2 "-" $3 "\t" $4}}' {input} > {output.id2gene} 2> {log}
		cut -f1 {output.id2gene} | tail -n +2 > {output.id}
		"""

rule createNormalPromoter:
	input: "output/{sample}/3_cancer/tfbs/id.txt"
	output: 
		temp("output/{sample}/3_cancer/tfbs/ref.promoter.fa")
	params:
		ref = genome_path
	singularity: "docker://quay.io/biocontainers/samtools:1.16.1--h6899075_1"
	log: "output/{sample}/log/createNormalPromoter.log"
	shell:
		"samtools faidx -r {input} {params.ref} > {output} 2> {log}"

rule mutatePromoter:
	input: 
		ref_promoter = "output/{sample}/3_cancer/tfbs/ref.promoter.fa",
		id = "output/{sample}/3_cancer/tfbs/id.txt",
		phase_vcf = lambda w: config["samples"][w.sample]["phase"]
	output: temp("output/{sample}/3_cancer/tfbs/allele{num}_promoter.fa")
	singularity: "docker://quay.io/biocontainers/bcftools:1.15--h0ea216a_2"
	log: "output/{sample}/log/mutatePromoter_{num}.log"
	shell:
		"""
		cat {input.ref_promoter} | bcftools consensus {input.phase_vcf} -H {wildcards.num}pIu > {output} 2> {log}
		"""

rule scanMotif:
	input: "output/{sample}/3_cancer/tfbs/allele{num}_promoter.fa"
	output: "output/{sample}/3_cancer/tfbs/allele{num}_motif.txt"
	params:
		motifFile="annotation/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
	singularity: "docker://quay.io/biocontainers/meme:5.4.1--py310pl5321hb021246_2"
	log: "output/{sample}/log/scanMotif_{num}.log"
	shell:
		"""
		fimo --text {params.motifFile} {input} > {output} 2> {log}
		"""
rule compareTFBS:
	input: 
		a1 = "output/{sample}/3_cancer/tfbs/allele1_motif.txt",
		a2 = "output/{sample}/3_cancer/tfbs/allele2_motif.txt",
		expression_matrix = getExpressionMatrix,
		id2gene = "output/{sample}/3_cancer/tfbs/id2gene.txt"
	output: "output/{sample}/3_cancer/tfbs/motifDiff.tsv"
	singularity: "docker://glenn032787/ase_rcontainer:2.0"
	params: 
		tf_list = "annotation/human_mono_motifs.tsv"
	log: "output/{sample}/log/compareTFBS.log"
	shell:
		"""
		Rscript scripts/compareTFBS.R \
			--allele1={input.a1} \
			--allele2={input.a2} \
			--expression_matrix={input.expression_matrix} \
			--id2gene={input.id2gene} \
			--tf={params.tf_list} \
			--min=10 \
			--outdir=output/{wildcards.sample}/3_cancer/tfbs \
			--sample={wildcards.sample} 2> {log}
		"""

### -------------------------------------------------------------------
### Cancer analysis - Stop gain/loss mutation
### -------------------------------------------------------------------

rule annotatePhase:
	input: lambda w: config["samples"][w.sample]["phase"]
	output: temp("output/{sample}/3_cancer/stopVar/phase.annotate.vcf")
	singularity: "docker://quay.io/biocontainers/snpeff:5.0--hdfd78af_1"
	params:
		genome = config["annotationPath"]["snpEff_genomeName"],
		snpEff_config = config["annotationPath"]["snpEff_config"],
		snpEff_datadir = config["annotationPath"]["snpEff_datadir"]
	log: "output/{sample}/log/annotatePhase.log"
	shell:
		"""
		mkdir -p output/{wildcards.sample}/3_cancer/stopVar

		snpEff -Xmx16g \
			-v {params.genome} \
			-c {params.snpEff_config} \
			-dataDir {params.snpEff_datadir} \
			-noStats \
			{input} > {output} #2> {log}

		"""	

rule oneLine:
	input: "output/{sample}/3_cancer/stopVar/phase.annotate.vcf"
	output: temp("output/{sample}/3_cancer/stopVar/phase.oneline.annotate.vcf")
	log: "output/{sample}/log/oneLine.log"
	shell:
		"""
		cat {input} | scripts/vcfEffOnePerLine.pl > {output} 2> {log}
		"""

rule snpSiftPhase:
	input: "output/{sample}/3_cancer/stopVar/phase.oneline.annotate.vcf"
	output: "output/{sample}/3_cancer/stopVar/phase.annotate.tsv"
	singularity: "docker://quay.io/biocontainers/snpsift:5.1d--hdfd78af_0"
	log: "output/{sample}/log/snpSiftPhase.log"
	shell:
		"""
		SnpSift -Xmx16g \
			extractFields {input} \
			ANN[0].GENE ANN[0].EFFECT GEN[0].GT FILTER > {output} 2> {log}
		"""

rule getStopMutation:
	input: "output/{sample}/3_cancer/stopVar/phase.annotate.tsv"
	output: "output/{sample}/3_cancer/stopVar/stop_variant.tsv"
	singularity: "docker://glenn032787/ase_rcontainer:2.0"
	log: "output/{sample}/log/getStopMutation.log"
	shell:
		"""
		Rscript scripts/stopVariant.R \
			--annotation={input} \
			--outdir=output/{wildcards.sample}/3_cancer/stopVar &> {log}
		"""	


### -------------------------------------------------------------------
### Cancer analysis - Somatic Mutation
### -------------------------------------------------------------------

rule getGeneSlop:
	input: 
		gene = "output/{sample}/3_cancer/raw/gene_annotation.bed",
		length = "output/{sample}/3_cancer/raw/genome.length"
	output: "output/{sample}/3_cancer/somatic/geneSlop.bed"
	singularity: "docker://quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
	log: "output/{sample}/log/getGeneSlop.log"
	shell:
		"""
		mkdir -p output/{wildcards.sample}/3_cancer/somatic
		bedtools slop -l 5000 -r 1000 -i {input.gene} -g {input.length} > {output} 2> {log}
		"""


rule somaticIntersect:
	input:
		genes = "output/{sample}/3_cancer/somatic/geneSlop.bed",
		variants = lambda w: config["samples"][w.sample][w.variant]
	output: "output/{sample}/3_cancer/somatic/{variant}.bed"
	singularity: "docker://quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6"
	log: "output/{sample}/log/somatic_{variant}_intersect.log"
	shell:
		"""
		bedtools intersect -a {input.genes} -b {input.variants} > {output} 2> {log}
		"""


### -------------------------------------------------------------------
### Cancer analysis - Summary Table
### -------------------------------------------------------------------

# Functions to check input

def checkCNV(wildcards):
	if config['cancer_analysis'] and "cnv" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["cnv"] != None:
		return "output/{sample}/3_cancer/intersect/cnv_intersect.bed"
	else:
		return []

def checkMethyl(wildcards):
	if config['cancer_analysis'] and "methyl" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["methyl"] != None:
		return "output/{sample}/3_cancer/intersect/methyl_intersect.bed"
	else:
		return []

def checkTFBS(wildcards):
	if config['cancer_analysis'] and "phase" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["phase"] != None:
		return "output/{sample}/3_cancer/tfbs/motifDiff.tsv"
	else:
		return []

def checkStopVar(wildcards): 
	if config['cancer_analysis'] and "phase" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["phase"] != None:
		return "output/{sample}/3_cancer/stopVar/stop_variant.tsv"
	else:
		return []

def checkSomaticSnv(wildcards):
	if config['cancer_analysis'] and "somatic_snv" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["somatic_snv"] != None:
		return "output/{sample}/3_cancer/somatic/somatic_snv.bed"
	else:
		return []	

def checkSomaticIndel(wildcards):
	if config['cancer_analysis'] and "somatic_indel" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["somatic_indel"] != None:
		return "output/{sample}/3_cancer/somatic/somatic_indel.bed"
	else:
		return []

def checkCancerAnalysis(wildcards):
	if config['cancer_analysis']:
		return "annotation/cancer_gene.txt"
	else: 
		return []

def checkTissue(wildcards):
	if config['cancer_analysis'] and "tissue" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["tissue"] != None and config["samples"][wildcards.sample]["tissue"] in possibleTissue:
		return config["samples"][wildcards.sample]["tissue"]
	else:
		return []
	

rule summaryTableCancer:
	input:
		cnv = checkCNV,
		methyl = checkMethyl, 
		tfbs = checkTFBS,
		stopVar = checkStopVar,
		snv = checkSomaticSnv,
		indel = checkSomaticIndel,
		ase = "output/{sample}/2_mBASED/MBASED_expr_gene_results.txt" 
	output: "output/{sample}/summaryTable.tsv"
	singularity: "docker://glenn032787/ase_rcontainer:2.0"
	params:
		cancer=checkCancerAnalysis,
		normal="annotation/phaserNormalASE.tsv",
		tissue=checkTissue
	log: "output/{sample}/log/summaryTable.log"
	shell:
		"""
		Rscript scripts/summaryTable.R \
			--cnv={input.cnv} \
			--methyl={input.methyl} \
			--tfbs={input.tfbs} \
			--stop={input.stopVar} \
			--snv={input.snv} \
			--indel={input.indel} \
			--ase={input.ase} \
			--sample={wildcards.sample} \
			--cancer={params.cancer} \
			--tissue={params.tissue} \
			--normal={params.normal} \
			--outdir=output/{wildcards.sample} 2> {log}
		"""


### -------------------------------------------------------------------
### Cancer analysis - Figures
### -------------------------------------------------------------------

rule cancerFigures:
	input: "output/{sample}/summaryTable.tsv"
	output: "output/{sample}/figures/aseCause.pdf"
	singularity: "docker://glenn032787/ase_rcontainer:2.0"
	log: "output/{sample}/log/cancerFigures.log" 
	shell:
		"""
		mkdir -p output/{wildcards.sample}/figures/tables
	
		Rscript scripts/cancerFigures.R \
			--summary={input} \
			--outdir=output/{wildcards.sample}/figures \
			--sample={wildcards.sample} &> {log}
		"""

### -------------------------------------------------------------------
### Karyogram Figure
### -------------------------------------------------------------------

def checkCNV_karyogram(wildcards):
	if config['cancer_analysis'] and "cnv" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["cnv"] != None:
		return config["samples"][wildcards.sample]["cnv"]
	else:
		return []

def checkMethyl(wildcards):
	if config['cancer_analysis'] and "methyl" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["methyl"] != None:
		return config["samples"][wildcards.sample]["methyl"]
	else:
		return []

rule karyogram:
	input:
		cnv = checkCNV_karyogram,
		dmr = checkMethyl,
		ase = "output/{sample}/summaryTable.tsv",
		centromere = config["centromere"][config["genome_name"]],
		chromSize = "output/{sample}/3_cancer/raw/genome.length",
		annotation = "annotation/biomart_ensembl100_GRCh38.sorted.bed"	
	output:
		"output/{sample}/figures/karyogram.pdf"
	singularity: "docker://glenn032787/ase_rcontainer:2.0"
	shell:
		"""
		Rscript scripts/karyogramFigure.R \
        		--chromSize={input.chromSize} \
        		--centPos={input.centromere} \
        		--cna={input.cnv} \
        		--dmr={input.dmr} \
        		--ase={input.ase} \
        		--genes={input.annotation} \
        		--out=output/{wildcards.sample}/figures/karyogram
		"""
	




	
