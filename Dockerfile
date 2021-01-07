# The format of the dockerfile can be found here: https://docs.docker.com/engine/reference/builder/
# See https://hub.docker.com/r/rocker/r-ver
FROM rocker/r-ver:4.0.3
MAINTAINER Bret Barnes "bbarnes@illumina.com"

# This is required by some sesame dependency
RUN apt-get update
RUN apt-get install zlib1g-dev libxml2 vim

# To install R packages use the following, I just pulled these real fast but you might need more
RUN install2.r dbplyr optparse tidyverse purrr tidyr plyr stringr readr glue ggplot2 tibble matrixStats scales doParallel R.utils scales GGally RcppArmadillo 

# You can't install bioconductor packages directly cause they decided on some stupid repo format - this is a fix around that
RUN install2.r littler BiocManager

RUN Rscript /usr/local/lib/R/site-library/littler/examples/installBioc.r Biostrings
RUN Rscript /usr/local/lib/R/site-library/littler/examples/installBioc.r data.table
RUN Rscript /usr/local/lib/R/site-library/littler/examples/installBioc.r AnnotationHub
RUN Rscript /usr/local/lib/R/site-library/littler/examples/installBioc.r ExperimentHub

RUN Rscript /usr/local/lib/R/site-library/littler/examples/installBioc.r sesame

COPY . /opt/Infinium_Methylation_Workhorse
WORKDIR /opt/Infinium_Methylation_Workhorse

# Note the actual R script will be in /opt/Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R
# Args should be passed in. Might need to rework your .sh wrapper (or just ditch it)
CMD [ Rscript /opt/Infinium_Methylation_Workhorse/scripts/R/swifthoof/swifthoof_main.R ]

# End of file
