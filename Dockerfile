FROM nfcore/base:1.10.2
LABEL authors="Barry Digby" \
      description="Docker image containing tools for circRNA analysis"
      
COPY environment.yml /

RUN conda env create -f environment.yml && conda clean -a

ENV PATH /opt/conda/envs/circrna/bin:$PATH

RUN sed -i '126s/^/#/' /opt/conda/envs/circrna/lib/python2.7/site-packages/CIRIquant/main.py
