FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y --no-install-recommends apt-transport-https software-properties-common dirmngr gpg-agent
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

RUN apt-get update && apt-get install -y --no-install-recommends build-essential gfortran libblas-dev liblapack-dev libz-dev bwa wget r-base libpng-dev imagemagick

RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('fitdistrplus')"
RUN Rscript -e "install.packages('classInt')"
RUN Rscript -e "install.packages('actuar')"
RUN Rscript -e "install.packages('stringr')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('PRROC')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('ComplexHeatmap')"

RUN mkdir /root/data
ADD data /root/data/
RUN mkdir /root/doc
RUN mkdir /root/src
ADD src /root/src/

RUN g++ -Wall -std=c++17 -I /root/src/alglib_include -O3 -DDEBUG -I /root -lz -pthread /root/src/Projects/Repare/ccaGatherQC.cpp -o /usr/local/bin/ccaGatherQC
RUN g++ -Wall -std=c++17 -I /root/src/alglib_include -O3 -DDEBUG -I /root -lz -pthread /root/src/Projects/Repare/cleanUpCca.cpp -o /usr/local/bin/cleanUpCca
RUN g++ -Wall -std=c++17 -I /root/src/alglib_include -O3 -DDEBUG -I /root -lz -pthread /root/src/Projects/Repare/crisprCountsAnalysis.cpp -o /usr/local/bin/crisprCountsAnalysis
RUN g++ -Wall -std=c++17 -I /root/src/alglib_include -O3 -DDEBUG -I /root -lz -pthread /root/src/Projects/Repare/barplotGenesFromCcaFoldchange.cpp -o /usr/local/bin/barplotGenesFromCcaFoldchange
RUN g++ -Wall -std=c++17 -I /root/src/alglib_include -O3 -DDEBUG -I /root -lz -pthread /root/src/Projects/Repare/convertCrisprFastqsToReadCountFile.cpp -o /usr/local/bin/convertCrisprFastqsToReadCountFile
