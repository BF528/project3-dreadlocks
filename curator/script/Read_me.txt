#####Read_me#####

#working directory#
/projectnb2/bf528/users/dreadlocks/project_3/curator

#copy the files we need in our directory#
scp /project/bf528/project_3/toxgroups/toxgroup_6_rna_info.csv /projectnb2/bf528/users/dreadlocks/project_3/curator/sample

bash /projectnb2/bf528/users/dreadlocks/project_3/curator/script/get_sample.sh

#run fastqc on the fastq files
qsub /projectnb2/bf528/users/dreadlocks/project_3/curator/script/run_fastqc.qsub

#run STAR
qsub /projectnb2/bf528/users/dreadlocks/project_3/curator/script/run_star.qsub

#run multiqc
qsub /projectnb2/bf528/users/dreadlocks/project_3/curator/script/run_multiqc.qsub