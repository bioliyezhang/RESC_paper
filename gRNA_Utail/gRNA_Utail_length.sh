#!/bin/bash
for i in $(cat gRNA.list)
do
samtools view -H example.bam > header.sam
samtools view example.bam | grep -w $i > example_g-.sam
cat header.sam example_g-.sam > example_g.sam
samtools view -bS example_g.sam > example_g.bam
python gRNA_Utail_length.py example_g.bam
done
