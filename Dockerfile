FROM continuumio/miniconda3

WORKDIR /app

# Make RUN commands use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]


# Create the environment:
RUN git clone https://github.com/kblin/covid-spike-classification.git covid-spike-classification
RUN ls && pwd
#COPY covid-spike-classification /app
#WORKDIR covid-spike-classification
RUN conda env create -n csc -f covid-spike-classification/environment.yml
RUN conda activate csc
RUN pip install covid-spike-classification



# Initialize conda in bash config files:
RUN conda init bash

# Activate the environment, and make sure it's activated:
RUN echo "conda activate csc" > ~/.bashrc




# The code to run when container is started:
