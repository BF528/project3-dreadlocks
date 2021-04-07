#!/bin/bash

#$ -P bf528
#$ -cwd
#$ -j y
#$ -pe mpi_16_tasks_per_node 16

SAMPLES=('SRR1177997' 'SRR1177999' 'SRR1178002' 'SRR1178014' 'SRR1178021' 'SRR1178047' 'SRR1177963' 'SRR1177964' 'SRR1177965')
SAMPLE_FILE='/projectnb2/bf528/users/dreadlocks/project_3/curator/sample'

echo "Running job $JOB_ID"
echo "Started: $(date +%F)"
echo "Running in directory: $PWD"

for i in ${SAMPLES[@]}; 
do
    file1=('/project/bf528/project_3/samples/'$i'_1.fastq.gz')
    file2=('/project/bf528/project_3/samples/'$i'_2.fastq.gz')
    scp $file1 $file2 $SAMPLE_FILE
done

