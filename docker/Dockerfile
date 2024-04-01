FROM rocker/verse:4.3.3

RUN apt-get update -qq
RUN apt-get -y install cmake
RUN apt-get -y --no-install-recommends install libglu1-mesa-dev
RUN apt-get -y install libxml2-dev
RUN apt-get -y install ghostscript
RUN apt-get -y install libxt-dev
RUN apt-get -y install git
RUN apt-get -y install -y libpython3-dev
RUN apt-get -y install python3-pip

RUN pip3 install https://github.com/wichmann-lab/python-psignifit/zipball/master

RUN install2.r –error \
--deps TRUE \
rmarkdown \
reticulate \
knitr \
osfr \
bookdown \
rstan \
coda \
modeest \
pals \
pracma \
doParallel \
R.matlab \
bayestestR \
alphahull \
kableExtra
