FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get -y install curl wget python3 python3-pip r-mathlib liblapack3 g++ gfortran
RUN pip3 install pandas scipy

RUN wget http://wp.unil.ch/sgg/files/2017/08/quicktest-1.1_bgen_v1.2.tar_.gz && \
tar -xzvf quicktest-1.1_bgen_v1.2.tar_.gz

RUN cd quicktest-1.1_bgen_v1.2 && \
make clean && \
make -f Makefile

RUN apt-get update && apt-get install -y dstat atop

COPY format_quicktest_phenos.py format_quicktest_output.py /
