FROM bioboxes/biobox-minimal-base

ADD image /usr/local
RUN install.sh
ADD *.py /usr/local/bin/ 
ADD utils /usr/local/bin/utils
ADD requirements /requirements
RUN pip3 install -r /requirements/default.txt

ENV BIOBOX_EXEC /usr/local/bin/evaluate.sh
ENV TASKFILE /usr/local/share/Taskfile
ENV SCHEMA /usr/local/share/schema.yaml
