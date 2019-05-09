FROM bioconductor/release_core2:R3.6.0_Bioc3.9
MAINTAINER saimukund@wustl.edu

#add a few useful tools
RUN apt-get update -y && apt-get install -y --no-install-recommends \
    curl \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libbz2-dev \
    liblzma-dev \
    libpng-dev \
    bzip2

#add r packages
ADD rpackages.R /tmp/
RUN R -f /tmp/rpackages.R

#add R script
ADD CellMatch_Haemopedia.r /opt/
