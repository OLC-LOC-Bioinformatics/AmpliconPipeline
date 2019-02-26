# Dockerfile for AmpliconPipeline

# NOTE: This Dockerfile will successfully build into a usable container, though due to temporary files created by DADA2,
# the full pipeline cannot run to completion. DADA2 will by default drop temporary files into the Docker container and
# eventually run out of space. Until this is fixed, use the conda installation method instead.

FROM ubuntu:16.04

MAINTAINER Forest Dussault <forest.dussault@inspection.gc.ca>

ENV DEBIAN_FRONTEND noninteractive

# Install packages
RUN apt-get update -y -qq && apt-get install -y \
	python-dev \
	git \
	curl \
	wget \
	python3-pip \
	ttf-dejavu \
	nano

ENV PATH /usr/sbin:$PATH
RUN useradd -ms /bin/bash/ ubuntu
USER ubuntu

WORKDIR HOME
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/ubuntu/miniconda.sh
RUN bash /home/ubuntu/miniconda.sh -b -p /home/ubuntu/miniconda
ENV PATH /home/ubuntu/miniconda/bin:$PATH
RUN echo $PATH
RUN conda install -y python=3 \
	    && conda update conda

# Upgrade pip
RUN pip3 install --upgrade pip

# Install AmpliconPipeline
WORKDIR /home/ubuntu/
ENV PATH /home/ubuntu/AmpliconPipeline:$PATH
RUN git clone https://github.com/forestdussault/AmpliconPipeline.git
WORKDIR /home/ubuntu/AmpliconPipeline
RUN conda create --name AmpliconPipeline --file requirements.txt

# Retrieve the classifier
RUN mkdir /home/ubuntu/AmpliconPipeline/classifiers
RUN curl -L https://ndownloader.figshare.com/files/10970087 -o classifiers/99_V3V4_Silva_naive_bayes_classifier.qza

# Setup data directory
RUN mkdir /home/ubuntu/AmpliconPipeline/data

# Set the language to use utf-8 encoding
ENV LANG C.UTF-8

# Building the image
# docker build --no-cache -t "ampliconpipeline:v0.9.1" .

# Pushing the image to Dockerhub
# docker push forestdussault/ampliconpipeline:v0.9.1

# Running the image interactively (mounting volume)
# docker run -it --rm --mount source=ampliconpipeline-vol,target=/data forestdussault/ampliconpipeline:v0.9
# docker run -it --rm --mount src=`pwd` forestdussault/ampliconpipeline:v0.9

# Running the image interactively (binding volume):
# docker run -it --rm -v /mnt/nas:/mnt/nas ampliconpipeline:v0.9.1

# Example command
# python ampliconpipeline.py -i /mnt/nas/MiSeq_Amplicon/171206 -o /mnt/nas/Forest/ampliconpipeline_docker_test -m /mnt/nas/Forest/ampliconpipeline_docker_test.tsv -f -tlf 18 -tlr 10