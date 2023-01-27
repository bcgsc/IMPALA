# Written by Glenn Chang
# Date: January 27th, 2023
# Utility function to help control the logic of the snakemake workflow

################
# Check if RNA bam file or R1, R2 RNA reads is used as input
################

def getAlignment(wildcards):
	"""
	Check if RNA bam file is exist, if so return RNA bam file (method 2). If not, run rules to obtain RNA bam file from R1, R2 reads (method 1)
	"""
	if "rna" not in config["samples"][wildcards.sample] or config["samples"][wildcards.sample]["rna"] == None:
		return alignmentFile
	return config["samples"][wildcards.sample]["rna"]

def getExpressionMatrix(wildcards):
	"""
	Check if RNA bam is avalible, if not run rule to get expression matrix from R1 and R2 RNA reads..
	"""
	if "rna" not in config["samples"][wildcards.sample] or  config["samples"][wildcards.sample]["rna"] == None:
		return "output/{sample}/0_alignment/expression_matrix.tsv"
	return rpkm_path

################
# Check optional input for tumor content
################
def getTumorContent(wildcards):
	"""
	Return tumor content if avalible in config file, if not return 1.0 (assume 100% tumor content)
	"""
	if "tumorContent" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["tumorContent"] != None:
		return config["samples"][wildcards.sample]["tumorContent"]
	else:
		return 1.0

###################
# Check optional inputs for summaryTableCancer rule. Return empty list if not avalible
###################

def checkCNV(wildcards):
	"""
	Check and return CNV data if it is provided in config file and cancer_analysis is True, if not return empty list
	"""
	if config['cancer_analysis'] and "cnv" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["cnv"] != None:
		return "output/{sample}/3_cancer/intersect/cnv_intersect.bed"
	else:
		return []

def checkMethyl(wildcards):
	"""
        Check and return methylation data if it is provided in config file and cancer_analysis is True, if not return empty list
        """
	if config['cancer_analysis'] and "methyl" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["methyl"] != None:
		return "output/{sample}/3_cancer/intersect/methyl_intersect.bed"
	else:
		return []

def checkTFBS(wildcards):
	"""
	Check if phasing is provided in config file and cancer_analysis is True, if so run rules to obtain allelic TFBS, else return empty list
	"""	
	if config['cancer_analysis'] and "phase" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["phase"] != None:
		return "output/{sample}/3_cancer/tfbs/motifDiff.tsv"
	else:
		return []

def checkStopVar(wildcards):
	"""
        Check if phasing is provided in config file and cancer_analysis is True, if so run rules to obtain stop variant, else return empty list
        """
	if config['cancer_analysis'] and "phase" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["phase"] != None:
		return "output/{sample}/3_cancer/stopVar/stop_variant.tsv"
	else:
		return []

def checkSomaticSnv(wildcards):
	"""
        Check and return somatic SNV data if it is provided in config file and cancer_analysis is True, if not return empty list
        """
	if config['cancer_analysis'] and "somatic_snv" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["somatic_snv"] != None:
		return "output/{sample}/3_cancer/somatic/somatic_snv.bed"
	else:
		return []

def checkSomaticIndel(wildcards):
	"""
        Check and return somatic indel data if it is provided in config file and cancer_analysis is True, if not return empty list
        """
	if config['cancer_analysis'] and "somatic_indel" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["somatic_indel"] != None:
		return "output/{sample}/3_cancer/somatic/somatic_indel.bed"
	else:
		return []

def checkCancerAnalysis(wildcards):
	"""
	Return cancer gene list if cancer_analysis is True
	"""
	if config['cancer_analysis']:
		return "annotation/cancer_gene.txt"
	else:
		return []

def checkTissue(wildcards):
	"""
	Return tissue type if avalible and if cancer_analysis is True, else return empty list
	"""
	if config['cancer_analysis'] and "tissue" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["tissue"] != None and config["samples"][wildcards.sample]["tissue"] in possibleTissue:
		return config["samples"][wildcards.sample]["tissue"]
	else:
		return []


######################
# Check optional inputs for karyogram rule. Return empty list if not avalible
######################

def checkCNV_karyogram(wildcards):
	"""
	Check and return CNV data if it is provided in config file and cancer_analysis is True, if not return empty list
	"""
	if config['cancer_analysis'] and "cnv" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["cnv"] != None:
		return config["samples"][wildcards.sample]["cnv"]
	else:
		return []

def checkMethyl(wildcards):
	"""
        Check and return methylation data if it is provided in config file and cancer_analysis is True, if not return empty list
        """
	if config['cancer_analysis'] and "methyl" in config["samples"][wildcards.sample] and config["samples"][wildcards.sample]["methyl"] != None:
		return config["samples"][wildcards.sample]["methyl"]
	else:
		return []
