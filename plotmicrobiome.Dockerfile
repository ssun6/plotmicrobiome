# Deployment path:  $DOCKER_DIR/plotmicrobiome.Dockerfile
FROM r-base:4.2.1


#install git
RUN mv /bin /bin.old && ln -s /usr/bin /bin && \
    mv /sbin /sbin.old && ln -s /usr/sbin /sbin && \
    mv /lib /lib.old && ln -s /usr/lib /lib && \
    mv /lib64 /lib64.old && ln -s /usr/lib64 /lib64 || true
#RUN apt-get -y install git
RUN cd /usr/local/bin && git clone https://github.com/ssun6/plotmicrobiome.git

#install R packages
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install("ggtree")"
RUN R -e "BiocManager::install("rhdf5")"
RUN R -e "BiocManager::install("shiny")"
RUN R -e "BiocManager::install("shinyjs")"