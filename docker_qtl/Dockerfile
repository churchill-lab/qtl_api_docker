FROM rocker/r-ver:3.3.2
MAINTAINER Matthew Vincent <mvincent@jax.org>

RUN apt-get update -qq && apt-get install -y \
  apt-utils \
  git-core \
  libssl-dev \
  libcurl4-gnutls-dev 

RUN R -e 'install.packages(c("devtools", "pryr", "yaml", "jsonlite", "dplyr", "data.table", "RcppEigen", "plumber", "RSQLite"))'
RUN R -e 'devtools::install_github(paste0("rqtl/qtl2", c("geno", "scan", "plot", "convert")))'
RUN R -e 'devtools::install_github("simecek/intermediate")'

ADD run_plumber.R /api/run_plumber.R
ADD start_server.sh /api/start_server.sh

EXPOSE 8000

CMD ["/api/start_server.sh"]

