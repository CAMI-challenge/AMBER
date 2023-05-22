FROM python:3.11.3-slim

ADD *.py /usr/local/bin/
ADD cami_amber /usr/local/bin/cami_amber
ADD cami_amber/utils /usr/local/bin/cami_amber/utils
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt
