FROM python:3.11.3-slim

ADD image /usr/local
ADD *.py /usr/local/bin/
ADD src /usr/local/bin/src
ADD src/utils /usr/local/bin/src/utils
ADD requirements.txt requirements.txt
RUN pip install -r requirements.txt
