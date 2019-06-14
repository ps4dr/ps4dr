FROM r-base

#FROM scipy
#RUN apt-get -y install r-base
#RUN pip install rpy2
#RUN apt-get -y install libcurl4-openssl-dev
#setup R configs
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN apt-get update \
	&& apt-get install -t unstable -y --no-install-recommends \
	    apt-utils \
	    openssl \
	    libcurl4-openssl-dev \
	    libxml2-dev

# Install development tools
RUN R -e 'install.packages("devtools")'

# Install this code from GitHub
RUN R -e 'library(devtools); install_github("ps4dr/ps4dr")'
