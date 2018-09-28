# This is an alpine linux installation

FROM alpine
MAINTAINER Simon Mudd (simon.m.mudd@ed.ac.uk) and Fiona Clubb (clubb@uni-potsdam.de)

# install essential packages
RUN apk upgrade -U && \
    apk update && \
    apk add --virtual build-dependencies build-base gcc wget git bash && \
    apk add --repository http://dl-cdn.alpinelinux.org/alpine/edge/main libressl2.7-libcrypto && \
    apk add gdal --update-cache --repository http://dl-cdn.alpinelinux.org/alpine/edge/testing && \
    apk add fftw-dev cmake && \
    rm -rf /var/cache/apk/* && \
    rm -rf /tmp/*

# update to avoid weird apk error 
RUN apk update

# Add an LSDTopoTools working directory
WORKDIR /LSDTopoTools/
