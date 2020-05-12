FROM python:3.7-slim

ADD image /usr/local
ADD *.py /usr/local/bin/
ADD src /usr/local/bin/src
ADD src/utils /usr/local/bin/src/utils
ADD requirements /requirements
RUN pip install -r /requirements/default.txt

ENV BIOBOX_EXEC /usr/local/bin/evaluate.sh
ENV TASKFILE /usr/local/share/Taskfile
ENV SCHEMA /usr/local/share/schema.yaml
