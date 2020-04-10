FROM ubuntu:18.04

MAINTAINER James Gallant <jgallant@sun.ac.za>

#main
RUN apt-get -qqy update -y
RUN apt-get -qqy install openjdk-8-jre -y
RUN apt-get -qqy install python3.8 -y 
RUN apt-get -qqy install python2.7 -y
RUN apt-get -qqy install curl -y
RUN apt-get -qqy install wget -y 
RUN apt-get -qqy install git -y
RUN apt-get -qqy install unzip -y
RUN apt-get -qqy install autoconf -y

RUN apt-get -qqy install libboost-all-dev
RUN apt-get -qqy install gcc -y
RUN apt-get -qqy install make -y 
RUN apt-get -qqy install libbz2-dev -y
RUN apt-get -qqy install zlib1g-dev -y 
RUN apt-get -qqy install libncurses5-dev -y
RUN apt-get -qqy install libncursesw5-dev -y
RUN apt-get -qqy install liblzma-dev -y
RUN apt-get -qqy install cython -y 
RUN apt-get -qqy install cython3 -y
RUN apt-get -qqy install python3-pip -y 
RUN apt-get -qqy install python-pip -y
RUN apt-get -qqy install libncursesw5-dev -y
RUN pip install --upgrade pip
RUN pip install pysam
RUN pip install numpy
RUN pip3 install pandas
RUN pip3 install pysam
RUN apt-get -qqy install samtools -y
RUN apt-get -qqy install bedtools -y
RUN apt-get -qqy install libcurl4-openssl-dev -y 
RUN apt-get -qqy install vcftools -y
RUN apt-get -qqy install libssl-dev -y 
RUN apt-get -qqy install bwa -y
RUN apt-get -qqy install bcftools -y
RUN apt-get -qqy install abacas -y
RUN apt-get -qqy install fastqc -y 
RUN apt-get -qqy install tabix -y
RUN apt-get -qqy install sambamba -y 
RUN apt-get -qqy install g++ -y
RUN apt-get -qqy install cmake -y
RUN apt-get -qqy install gawk -y 
RUN apt-get -qqy install gnuplot-nox -y
RUN apt-get -qqy install soapdenovo2 -y

WORKDIR /home/pegasus

COPY . .

RUN mkdir programs && \
    mv databases/1.jar ./programs/GenomeAnalsysTK.jar && \
    mv databases/2.tar.gz ./programs/novocraft.tar.gz && \
    cd ./programs && \
    tar -xvf novocraft.tar.gz

RUN wget -P ./programs/ "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip" && \
    wget -P ./programs/ "https://github.com/broadinstitute/picard/releases/download/2.12.2/picard.jar"  && \
    unzip ./programs/Trimmomatic-0.36.zip  && \
    mv Trimmomatic-0.36 ./programs  && \
    rm -r ./programs/Trimmomatic-0.36.zip

RUN mkdir programs/samblaster

RUN git clone --recursive "https://github.com/GregoryFaust/samblaster.git" programs/samblaster && \
    cd programs/samblaster/ && \
    make

RUN export PATH=$PATH:./programs/samblaster

RUN cd programs && \
    git clone --recursive https://github.com/hall-lab/svtyper.git -b v0.1.4

RUN mkdir programs/lumpy-sv && \
    git clone --recursive "https://github.com/arq5x/lumpy-sv.git" programs/lumpy-sv/ && \
    cd programs/lumpy-sv && \
    make 

RUN mkdir programs/delly && \
    git clone --recursive "https://github.com/dellytools/delly.git" programs/delly && \
    cd programs/delly && \
    make all

RUN mkdir programs/svprops && \
    git clone --recursive "https://github.com/dellytools/svprops.git" programs/svprops && \
    cd programs/svprops && \
    make
	    





