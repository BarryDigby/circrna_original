FROM nfcore/base:1.10.2
LABEL authors="Barry Digby" \
      description="Docker image containing tools for circRNA analysis"

# install main packages:
RUN apt-get update; apt-get clean all;

RUN apt-get install --yes build-essential \
                          gcc-multilib \
                          apt-utils \
			  curl \
                          perl \
			  zip \
                          expat \
                          libexpat-dev
			  
RUN apt-get install --yes cpanminus

RUN apt-get install --yes libxml-libxml-perl \
			  libxml-dom-xpath-perl \
			  libxml-libxml-simple-perl \
			  libxml-dom-perl

RUN cpanm CPAN::Meta Statistics::Lite Bio::TreeIO

# TargetScan Executables
RUN curl --output ./targetscan_70.zip http://www.targetscan.org/vert_72/vert_72_data_download/targetscan_70.zip
RUN unzip targetscan_70.zip
RUN curl --output ./targetscan_70_BL_PCT.zip http://www.targetscan.org/vert_72/vert_72_data_download/targetscan_70_BL_PCT.zip
RUN unzip targetscan_70_BL_PCT.zip
RUN curl --output ./TargetScan7_context_scores.zip http://www.targetscan.org/vert_72/vert_72_data_download/TargetScan7_context_scores.zip
RUN unzip TargetScan7_context_scores.zip 
RUN mv targetscan_70.pl /usr/bin/
RUN mv TargetScan7_BL_PCT/targetscan_70_BL_bins.pl /usr/bin/
RUN mv TargetScan7_BL_PCT/targetscan_70_BL_PCT.pl /usr/bin/
RUN chmod 777 TargetScan7_context_scores/targetscan_70_context_scores.pl && mv TargetScan7_context_scores/targetscan_70_context_scores.pl /usr/bin/
RUN mv TargetScan7_context_scores/targetscan_count_8mers.pl /usr/bin/

## ViennaRNA
WORKDIR /usr/bin
RUN wget --no-check-certificate https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.15.tar.gz
RUN tar -zxvf ViennaRNA-2.4.15.tar.gz
WORKDIR /usr/bin/ViennaRNA-2.4.15/
RUN ./configure --prefix=/usr/bin/ViennaRNA
RUN make install
RUN mv /usr/bin/ViennaRNA/bin/* /usr/bin

# Conda ENV    
WORKDIR /  
COPY environment.yml /
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/circrna/bin:$PATH
# comment CIRIquant line that attempts to run os.chmod() 
RUN sed -i '126s/^/#/' /opt/conda/envs/circrna/lib/python2.7/site-packages/CIRIquant/main.py

# DCC
WORKDIR /usr/src/app
RUN wget --no-check-certificate https://github.com/dieterich-lab/DCC/archive/v0.4.8.tar.gz
RUN tar -xvf v0.4.8.tar.gz
WORKDIR /usr/src/app/DCC-0.4.8
# remove --user or else scripts installed to /root/ 
RUN python setup.py install

# UROBORUS
WORKDIR /usr/src/app
RUN wget --no-check-certificate https://github.com/WGLab/UROBORUS/archive/v2.0.0.tar.gz
RUN tar -xvf v2.0.0.tar.gz
WORKDIR /usr/src/app/UROBORUS-2.0.0/bin
RUN chmod 777 UROBORUS.pl && cp UROBORUS.pl /opt/conda/envs/circrna/bin

# find_circ
WORKDIR /usr/src/app
RUN wget --no-check-certificate http://www.circbase.org/download/find_circ.tar.gz
RUN tar -xvf find_circ.tar.gz
RUN cp *.py /opt/conda/envs/circrna/bin

# circRNA_finder 
RUN wget --no-check-certificate https://github.com/orzechoj/circRNA_finder/archive/v1.2.tar.gz
RUN tar -xvf v1.2.tar.gz
WORKDIR /usr/src/app/circRNA_finder-1.2
RUN cp *.pl /opt/conda/envs/circrna/bin
RUN chmod 777 filterCirc.awk && cp filterCirc.awk /opt/conda/envs/circrna/bin

# R markdown to HTML
RUN R -e "install.packages('DT',dependencies=TRUE, repos='http://cran.rstudio.com/')"
