FROM umccr/umccrise:0.15.10

RUN R -e "install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')"
RUN R -e "BiocManager::install('TxDb.Hsapiens.UCSC.hg38.knownGene')"

