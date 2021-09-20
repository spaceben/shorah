FROM ubuntu:20.04

RUN apt-get update -y && \
    apt-get install -y autoconf automake m4 python3 python3-pip \
        libhts-dev libboost-math-dev libboost-random-dev zlib1g-dev && \ 
    DEBIAN_FRONTEND=noninteractive apt-get install -y pkg-config 

SHELL ["/bin/bash", "-c"]

COPY . .

RUN autoreconf -vif -I m4

RUN pip3 install Biopython numpy
RUN ./configure 

#ENTRYPOINT ["make", "-j1", "distcheck"]
ENTRYPOINT echo "$PWD"
