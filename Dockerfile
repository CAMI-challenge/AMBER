FROM bioboxes/biobox-minimal-base

RUN apt update
RUN apt install -y python python-pip python-dev 
ADD *.py /usr/local/bin/ 
ADD utils /usr/local/bin/utils
ADD requirements /requirements
RUN pip install -r /requirements/default.txt
ADD image /usr/local

ENV BIOBOX_EXEC /usr/local/bin/evaluate.sh
ENV TASKFILE /usr/local/share/Taskfile
ENV SCHEMA /usr/local/share/schema.yaml
