## Load config values
configfile: "config/defaults.yaml"
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"

# path to the reference genome fasta
genome_path = config["genome"][config["genome_name"]]

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
	    expand("output/{sample}/mBASED",sample=sample_ids)

### -------------------------------------------------------------------
### Call and filter the DNA SNVs
### -------------------------------------------------------------------

rule dna_snv_calling:
    input:
        bam = lambda w: config["samples"][w.sample]["dna"],
        ref = genome_path
    output:
        "output/{sample}/StrelkaDNA/results/variants/genome.S1.vcf.gz"
    conda: "config/strelka.yaml"
    threads: 20
    shell:
        """
        configureStrelkaGermlineWorkflow.py --bam={input.bam} --referenceFasta={input.ref} --rna --runDir=output/{wildcards.sample}/StrelkaDNA
        output/{wildcards.sample}/StrelkaDNA/runWorkflow.py -m local -j {threads}
        """

rule dna_snv_filt:
    input:
        vcf = "output/{sample}/StrelkaDNA/results/variants/genome.S1.vcf.gz"
    output:
        "output/{sample}/dna.het.pass.snps.vcf.gz"
    conda: "config/ase-env.yaml"
    shell:
        "zcat {input.vcf} | grep -E '(PASS|#)' | grep -E '(0/1|#)' | awk '/^#/||length($4)==1 && length($5)==1' | bgzip > {output}"

rule dna_snv_index:
    input:
        vcf="output/{sample}/dna.het.pass.snps.vcf.gz"
    output:
        "output/{sample}/dna.het.pass.snps.vcf.gz.tbi"
    conda: "config/ase-env.yaml"
    shell:
        "tabix {input.vcf}"

### -------------------------------------------------------------------
### Call and filter the RNA SNVs
### -------------------------------------------------------------------

rule rna_snv_calling:
    input:
        bam = lambda w: config["samples"][w.sample]["rna"],
        vcf = "output/{sample}/dna.het.pass.snps.vcf.gz",
        ref = genome_path
    output:
        "output/{sample}/StrelkaRNA/results/variants/genome.S1.vcf.gz"
    conda: "config/strelka.yaml"
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
    shell:
        "zcat {input.vcf} | grep -E '(PASS|#)' | bgzip > {output}"

rule rna_snv_index:
    input:
        vcf = "output/{sample}/rna.forceGT.pass.vcf.gz"
    output:
        "output/{sample}/rna.forceGT.pass.vcf.gz.tbi"
    conda: "config/ase-env.yaml"
    shell:
        "tabix {input.vcf}"

rule intersect:
    input:
        vcf1 = "output/{sample}/dna.het.pass.snps.vcf.gz",
        vcf2 = "output/{sample}/rna.forceGT.pass.vcf.gz"
    output:
        "output/{sample}/rna.isec.dna.snps.vcf"
    conda: "config/ase-env.yaml"
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
    shell:
        "bgzip {input}"

### -------------------------------------------------------------------
### Run MBASED
### -------------------------------------------------------------------

rule mbased:
    input:
        phase = lambda w: config["samples"][w.sample]["phase"],
        vcf = "output/{sample}/rna.isec.dna.snps.vcf.gz",
        rpkm = rpkm_path
    output:
        directory("output/{sample}/mBASED")
    shell:
        "scripts/mbased.R --phase={input.phase} --rna={input.vcf} --sample={wildcards.sample} --rpkm={input.rpkm} --min=1 --outdir={output}"