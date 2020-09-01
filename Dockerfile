FROM nfcore/base:1.10.2
LABEL authors="Barry Digby" \
      description="Docker image containing tools for circRNA analysis"
      
COPY environment.yml /

RUN conda env create -f environment.yml && conda clean -a

ENV PATH /opt/conda/envs/circrna/bin:$PATH
