FROM nfcore/base:1.10.2
LABEL authors="Barry Digby" \
      description="Docker image containing tools for circRNA analysis"

#Conda ENV      
COPY environment.yml /
RUN conda env create -f environment.yml && conda clean -a
ENV PATH /opt/conda/envs/circrna/bin:$PATH
# comment CIRIquant line that attempts to run os.chmod() 
RUN sed -i '126s/^/#/' /opt/conda/envs/circrna/lib/python2.7/site-packages/CIRIquant/main.py

#DCC
WORKDIR /usr/src/app
RUN wget --no-check-certificate https://github.com/dieterich-lab/DCC/archive/v0.4.8.tar.gz
RUN tar -xvf v0.4.8.tar.gz
WORKDIR /usr/src/app/DCC-0.4.8
# remove --user or else scripts installed to /root/ 
RUN python setup.py install

#UROBORUS
WORKDIR /usr/src/app
RUN wget --no-check-certificate https://github.com/WGLab/UROBORUS/archive/v2.0.0.tar.gz
RUN tar -xvf v2.0.0.tar.gz
WORKDIR /usr/src/app/UROBORUS-2.0.0/bin
RUN chmod 777 UROBORUS.pl && cp UROBORUS.pl /opt/conda/envs/circrna/bin

#find_circ
WORKDIR /usr/src/app
RUN wget --no-check-certificate http://www.circbase.org/download/find_circ.tar.gz
RUN tar -xvf find_circ.tar.gz
RUN cp *.py /opt/conda/envs/circrna/bin

#circRNA_finder 
RUN wget --no-check-certificate https://github.com/orzechoj/circRNA_finder/archive/v1.2.tar.gz
RUN tar -xvf v1.2.tar.gz
WORKDIR /usr/src/app/circRNA_finder-1.2
RUN cp *.pl /opt/conda/envs/circrna/bin
RUN chmod 777 filterCirc.awk && cp filterCirc.awk /opt/conda/envs/circrna/bin
