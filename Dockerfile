#Requires the GATK jar file, and the Platypus folder to be placed in the current directory

FROM ubuntu:16.04

RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install build-essential -y
RUN apt-get install git -y
RUN apt-get install wget -y
RUN apt-get install unzip -y
RUN apt-get install zlib1g-dev -y
RUN apt-get install libncurses5-dev -y
RUN apt-get install bzip2 -y
RUN apt-get install libbz2-dev -y
RUN apt-get install --reinstall make -y
RUN apt-get install cmake -y
RUN apt-get install gfortran -y
RUN apt-get install -y libatlas-base-dev
RUN apt-get install software-properties-common python-software-properties -y

# install python and PIP
RUN apt-get install python2.7-dev -y
#RUN apt-get install --upgrade cython -y
RUN apt-get install python-setuptools -y
RUN apt-get install python-pip -y
RUN pip install --upgrade pip setuptools wheel cython
RUN apt-get install python-numpy -y

RUN \
  echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
  add-apt-repository -y ppa:webupd8team/java && \
  apt-get update && \
  apt-get install -y oracle-java8-installer && \
  rm -rf /var/lib/apt/lists/* && \
  rm -rf /var/cache/oracle-jdk8-installer

ENV JAVA_HOME /usr/lib/jvm/java-8-oracle

RUN apt-get clean -y
WORKDIR /opt/

RUN git clone --recursive git://github.com/ekg/freebayes.git
WORKDIR /opt/freebayes
RUN make
RUN make install
WORKDIR /opt/freebayes/vcflib
RUN make
WORKDIR /opt/

# samtools
RUN apt-get install git -y
RUN git clone git://github.com/samtools/htslib.git
RUN git clone git://github.com/samtools/bcftools.git
RUN git clone git://github.com/samtools/samtools.git
RUN cd htslib && make && make install
WORKDIR /opt/
RUN cd samtools && make && make install
WORKDIR /opt/
RUN cd bcftools && touch doc/bcftools.1 && make && make install
WORKDIR /opt/

# numpy and such
RUN pip install --upgrade numpy
RUN pip install --upgrade pandas
RUN pip install -U git+https://github.com/pysam-developers/pysam
#RUN pip install pybedtools
#RUN pip install bx-python

# bwa
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /opt/bwa
RUN make
WORKDIR /opt

# hap.py
RUN git clone https://github.com/Illumina/hap.py.git /opt/hap.py-source
WORKDIR /opt/hap.py-source
RUN python install.py /opt/hap.py --no-tests
WORKDIR /opt
RUN rm -rf /opt/hap.py-source

# platypus
RUN git clone https://github.com/andyrimmer/Platypus.git /opt/platypus
WORKDIR /opt/platypus
RUN make
WORKDIR /opt

# gatk 3.5
RUN mkdir /opt/gatk-3.5-0
COPY GenomeAnalysisTK-3.5-0-g36282e4.tar.bz2 /opt/gatk-3.5-0
WORKDIR /opt/gatk-3.5-0
RUN tar xjf GenomeAnalysisTK-3.5-0-g36282e4.tar.bz2 && rm GenomeAnalysisTK-3.5-0-g36282e4.tar.bz2

# gatk 3.6
RUN mkdir /opt/gatk-3.6-0
COPY GenomeAnalysisTK-3.6.tar.bz2 /opt/gatk-3.6-0
WORKDIR /opt/gatk-3.6-0
RUN tar xjf GenomeAnalysisTK-3.6.tar.bz2 && rm GenomeAnalysisTK-3.6.tar.bz2
WORKDIR /opt

# vcfval
RUN wget https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.6.1/rtg-tools-3.6.1-linux-x64.zip

RUN unzip rtg-tools-3.6.1-linux-x64.zip && rm rtg-tools-3.6.1-linux-x64.zip

# vt
RUN git clone https://github.com/atks/vt.git
RUN cd vt && make
WORKDIR /opt

# vgraph
RUN pip install -U git+https://github.com/bioinformed/vgraph.git

#varscan
RUN wget 'http://heanet.dl.sourceforge.net/project/varscan/VarScan.v2.3.9.jar'

RUN apt-get update
RUN apt-get install python-scipy -y
RUN apt-get install python-matplotlib -y
WORKDIR /opt/
RUN git clone https://github.com/goranrakocevic/varcomp.git
WORKDIR /opt/varcomp
RUN pip install .
