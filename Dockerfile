FROM rocker/rstudio:4.1.2

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update --fix-missing && \
  apt-get install -y \
    libgit2-dev \
    openssh-client \
    liblzma-dev \
    libbz2-dev \
    libxml2-dev \
    libglpk-dev \
    libxt6 \
    && apt-get clean \
	  && rm -rf /var/lib/apt/lists/*

RUN echo "R_LIBS=/usr/local/lib/R/host-site-library:\${R_LIBS}" > /usr/local/lib/R/etc/Renviron.site \
  && echo "options(defaultPackages=c(getOption('defaultPackages'),'BiocManager'))" >> /usr/local/lib/R/etc/Rprofile.site \
  && wget -O Renviron.bioc http://bioconductor.org/checkResults/devel/bioc-LATEST/Renviron.bioc \
  && cat Renviron.bioc | grep -o '^[^#]*' | sed 's/export //g' >>/etc/environment \
  && cat Renviron.bioc >> /usr/local/lib/R/etc/Renviron.site \
  && echo BIOCONDUCTOR_VERSION=${BIOCONDUCTOR_VERSION} >> /usr/local/lib/R/etc/Renviron.site \
  && echo BIOCONDUCTOR_DOCKER_VERSION=${BIOCONDUCTOR_DOCKER_VERSION} >> /usr/local/lib/R/etc/Renviron.site \
  && echo 'LIBSBML_CFLAGS="-I/usr/include"' >> /usr/local/lib/R/etc/Renviron.site \
  && echo 'LIBSBML_LIBS="-lsbml"' >> /usr/local/lib/R/etc/Renviron.site \
  && rm -rf Renviron.bioc

ENV LIBSBML_CFLAGS="-I/usr/include"
ENV LIBSBML_LIBS="-lsbml"

COPY r_packages.R /home/rstudio/r_packages.R
RUN R -e "install.packages('BiocManager')" \
  && R -e "source('/home/rstudio/r_packages.R')"

USER rstudio
WORKDIR /home/rstudio
ENV PATH="/home/rstudio/conda/envs/reticulate/bin:/home/rstudio/conda/bin:$PATH"
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/rstudio/miniconda.sh \
  && /bin/bash /home/rstudio/miniconda.sh -b -p /home/rstudio/conda \
  && rm /home/rstudio/miniconda.sh \
  && /home/rstudio/conda/bin/conda init \
  && conda update -n base -c defaults conda -y \
  && conda install -c conda-forge mamba \
  && mamba init bash \
  && mamba create -n r-reticulate -c conda-forge python=3.8 leidenalg numpy python-igraph pandas \
  && mamba clean --all -y \
  && echo "export RETICULATE_AUTOCONFIGURE=0" >> /home/rstudio/.bashrc \
  && echo "export RETICULATE_AUTOCONFIGURE=0" >> /home/rstudio/.Renviron
ENV RETICULATE_PYTHON /home/rstudio/conda/envs/r-reticulate/bin

USER root
CMD ["/init"]