# Getting base image Ubuntu
FROM rocker/tidyverse:4.0.3

MAINTAINER Glenn <glchang@bcgsc.ca>

# Instal R CRAN package
RUN apt-get update && \
	apt-get install -y build-essential libglpk40 && \
	install2.r --error --skipinstalled \
	optparse RColorBrewer networkD3 htmlwidgets reshape2 ggrepel ggsci prob stats cowplot && \
	rm -rf /tmp/downloaded_packages

# Install bioconductor packages
RUN install2.r --error --skipinstalled BiocManager && \
	R -e 'BiocManager::install(ask = F)' && \
	R -e 'BiocManager::install(c("MBASED", "chromPlot", "SummarizedExperiment", ask = F))'
