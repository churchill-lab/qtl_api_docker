FROM rocker/hadleyverse:3.3.2
MAINTAINER Matthew Vincent <mvincent@jax.org>

RUN alias make='make -j'

RUN apt-get update -qq && apt-get install -y \
  apt-utils \
  git-core \
  libssl-dev \
  libcurl4-gnutls-dev 

RUN R -e 'install.packages(c("devtools", "yaml", "jsonlite", "dplyr", "data.table", "RcppEigen", "plumber", "RSQLite"))'
RUN R -e 'devtools::install_github(paste0("rqtl/qtl2", c("geno", "scan", "plot", "convert")))'
RUN R -e 'devtools::install_github("simecek/intermediate")'

ADD plumber_runner.R /opt/R/plumber_runner.R
ADD default.R /opt/R/default.R

EXPOSE 8000

ENTRYPOINT ["/usr/local/bin/Rscript", "/opt/R/plumber_runner.R"]
CMD ["/opt/R/default.R"]

