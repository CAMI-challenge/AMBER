#!/bin/bash
set -o errexit
set -o nounset
set -o xtrace

FASTA=$(biobox_args.sh 'select(has("fasta")) | map(.value) | join("")')

LABELS=$(biobox_args.sh 'select(has("labels")) | map(.value) | join("")')

PREDICTIONS=$(biobox_args.sh 'select(has("predictions")) | map(.value) | join("")')

TASK=$(fetch_task_from_taskfile.sh ${TASKFILE} $1)

OUTPUT_FILE="${OUTPUT}/metrics.txt"

eval $TASK

cat << EOF > ${OUTPUT}/biobox.yaml
version: 0.1.1
results:
  - name: Unsupervised Binning Evaluation
    type: tsv
    inline: false
    description: Unsupervised Binning Evaluation 
    value: metrics.txt
EOF
