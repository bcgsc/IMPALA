
## Allele specific expression
This workflow outputs the allele-specific expression using short-read RNA-seq and DNA-seq bams. The phasing information of the variants from tools such as WhatsApp can be provided to increase the performance of the tool (to-be an optional input but I don't know how to write that in snakemake - it is optional in the R script though)

### How to run
Clone the repository:
`git clone `

Right now conda environments are used, but it would be good to eventually change this to docker containers. This is the command to run:
`snakemake --use-conda -c 30`