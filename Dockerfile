# Requirements for all of my compbio containers.
#
FROM pytorch/pytorch AS compbio
RUN \
    apt-get update \
    && apt-get install -y \
        neovim \
        less \
        git \
        make \
        g++ \
        libz-dev \
        libtbb2 \
    && rm -rf /var/lib/apt/lists/*
RUN \
    conda install -y -c conda-forge \
        mamba \
    && conda clean -afy
RUN \
    mamba install -y -c conda-forge -c bioconda \
        cython \
        ipython \
        jupyter_contrib_nbextensions \
        jupyter_nbextensions_configurator \
        jupyterlab \
        lz4 \
        matplotlib \
        networkx \
        numpy \
        pandas \
        patsy \
        pigz \
        pip \
        python-graphviz \
        scikit-learn \
        scipy \
        seaborn \
        snakemake \
        statsmodels \
        tqdm \
        xarray \
        nodejs \
        rpy2 \
        xlrd openpyxl \
    && conda clean -afy
RUN \
    mamba install -y -c conda-forge \
        parallel \
    && (echo "will cite" | parallel --citation || true) \
    && conda clean -afy
RUN \
    mamba install -y -c conda-forge -c bioconda \
        biopython \
        bowtie2 \
        emboss \
        fastqc \
        fastuniq \
        hmmer \
        megahit \
        muscle \
        prodigal \
        prokka \
        pysam \
        samtools \
        scikit-bio \
        seqtk \
        diamond \
        bedtools \
        fasttree \
    && conda clean -afy
RUN \
    mamba install -y -c \
        conda-forge \
        pyro-ppl \
    && conda clean -afy
RUN \
    mamba install -y -c conda-forge \
        pymc3 \
        arviz \
    && conda clean -afy
RUN \
    mkdir -p src && cd src \
    && git clone https://github.com/zjshi/gt-pro gt-pro && cd gt-pro \
    && make && ln -s $PWD/GT_Pro /usr/bin/
RUN \
    mkdir -p src && cd src \
    && git clone https://github.com/lh3/gfa1 gfa1 && cd gfa1/misc \
    && make && ln -s $PWD/fastg2gfa /usr/bin/
RUN \
    mkdir -p src && cd src \
    && git clone https://github.com/najoshi/sickle sickle && cd sickle \
    && make && ln -s $PWD/sickle /usr/bin/
RUN \
    mkdir -p src && cd src \
    && git clone https://github.com/vsbuffalo/scythe scythe && cd scythe \
    && make all && ln -s $PWD/scythe /usr/bin/
RUN \
    mamba install -y -c conda-forge -c bioconda \
        trimmomatic \
    && conda clean -afy
RUN \
    mamba install -y -c conda-forge -c bioconda \
       r-vegan \
       r-dplyr \
       r-geepack \
       r-car \
    && conda clean -afy
RUN \
    mamba install -y -c conda-forge \
        unzip \
    && conda clean -afy
RUN \
    mamba install -y -c conda-forge \
        matplotlib-venn \
    && conda clean -afy
