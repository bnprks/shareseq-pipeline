# This base is ubunto 22.04 with fast R package installation 
# (downloads pre-compiled packages)
FROM eddelbuettel/r2u:22.04

# OS-package install
RUN apt-get update \
    && apt-get install -y parallel;

# Conda install
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh \
    && bash ~/miniconda.sh -b -p $HOME/miniconda \
    && ln -s $HOME/miniconda/bin/conda /usr/local/bin \
    && conda init \
    && conda install -y -c bioconda -c conda-forge \
        bedtools \
        bowtie2 \ 
        fastp \
        pip \
        rseqc \
        samtools \
        star \
        subread \
    ;

# Install umi_tools via pip
RUN $HOME/miniconda/bin/pip install umi_tools

# CRAN R packages (install.r is a shortcut from the base image)
RUN install.r \
        BiocManager \
        data.table \
        dplyr \
        Matrix \
        matrixStats \
        tibble \
        tidyr \
    && /usr/local/lib/R/site-library/littler/examples/installBioc.r \
        rhdf5 \
        DropletUtils \
    ;