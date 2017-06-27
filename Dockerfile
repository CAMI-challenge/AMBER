FROM ubuntu:16.04

RUN apt update
RUN apt install -y python python-pip 
ADD *.py /usr/local/bin/ 
ADD utils /usr/local/bin/utils
ADD requirements /requirements
RUN pip install -r /requirements/default.txt
