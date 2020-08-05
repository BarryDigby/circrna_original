FROM nfcore/base
LABEL authors="Barry Digby" \
      description="Docker image containing all requirements for circRNA analysis"
      
COPY environment.yml /
RUN conda env create -f /environment.yml python=2.7.15 && conda clean -a
ENV PATH /opt/conda/envs/circrna/bin:$PATH

COPY requirements.txt /
RUN pip install -r requirements.txt
