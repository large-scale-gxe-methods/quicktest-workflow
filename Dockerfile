FROM ubuntu

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install curl wget r-mathlib liblapack3 g++ gfortran

RUN wget http://wp.unil.ch/sgg/files/2017/08/quicktest-1.1_bgen_v1.2.tar_.gz && \
tar -xzvf quicktest-1.1_bgen_v1.2.tar_.gz

RUN cd quicktest-1.1_bgen_v1.2 && \
make clean && \
make -f Makefile




