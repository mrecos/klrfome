# get the base image, the rocker/verse has R, RStudio and pandoc
FROM rocker/verse

# required
MAINTAINER Your Name <mr.ecos@gmail.com>

COPY . /klrfome

# go into the repo directory
RUN . /etc/environment \

  # Install linux depedendencies here
  # e.g. need this for ggforce::geom_sina
  && sudo apt-get update \
  && sudo apt-get install libudunits2-dev -y \
  && sudo apt-get install libgeos-dev -y\

  # build this compendium package
  && R -e "install.packages('corrplot')"\
  && R -e "install.packages('dplyr')"\
  && R -e "install.packages('rdist')"\
  && R -e "install.packages('boot')"\
  && R -e "install.packages('pROC')"\
  && R -e "install.packages('raster')"\
  && R -e "install.packages('NLMR')"\
  && R -e "install.packages('rgeos')"\
  && R -e "install.packages('doParallel')"\
  && R -e "devtools::install('/klrfome', dep=TRUE)" \

 # render the manuscript into a docx, you'll need to edit this if you've
 # customised the location and name of your main Rmd file
  && R -e "rmarkdown::render('/klrfome/analysis/paper/paper.Rmd')"
