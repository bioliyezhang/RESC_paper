#!/bin/bash!

strain=$1
In_path=$2
Ref_gene=$3
UMI=$4
bin=$5
work_directory=$6
RNA_reverse=$7

mkdir -p ${work_directory}/0-split_fastq
mkdir ${work_directory}/1-1st_align
mkdir -p ${work_directory}/2-fastq_extract/split
mkdir ${work_directory}/3-2nd_align
mkdir -p ${work_directory}/4-UMI_merge/fastq/fail
mkdir ${work_directory}/4-UMI_merge/taf
mkdir -p ${work_directory}/5-select_final/taf
mkdir ${work_directory}/5-select_final/fastq
mkdir ${work_directory}/5-select_final/sam

## 0-split_fastq
split -l 2500000 ${In_path}/${strain}.fastq ${work_directory}/0-split_fastq/${strain}.fastq.
for gene in $(ls ${Ref_gene})
do
	for file in $(ls ${work_directory}/0-split_fastq/${strain}.fastq.*)
	do
		file_name=`echo $file | awk -F'split_fastq/' '{print $2}'`

		## 1. 1st align
		alignlib --in_ref ${Ref_gene}/${gene}/${gene}_maxicircle.fa 		\
	 	--in_lib ${file} 						\
	 	--out_prefix ${work_directory}/1-1st_align/${file_name}.${gene}.1st_Align.  	

		## 2. fastq extract
		python ${bin}/main/stranded_extract.py \
			-i ${work_directory}/1-1st_align/${file_name}.${gene}.1st_Align.mapped_reads.taf 		\
			-r ${file}							\
			-o ${work_directory}/2-fastq_extract/${strain}.${gene}		\
			-R ${RNA_reverse}
	done

	for strand in fwd rev
	do

		## 3. split
		split -l 2500000 ${work_directory}/2-fastq_extract/${strain}.${gene}.${strand}.fastq ${work_directory}/2-fastq_extract/split/${strain}.${gene}.${strand}.fastq.

		## 4. 2nd align
		for file in $(ls ${work_directory}/2-fastq_extract/split/${strain}.${gene}.${strand}.fastq.*)
		do
			file_name=`echo $file | awk -F'split/' '{print $2}'`
			alignlib --in_ref ${Ref_gene}/${gene}/${gene}_maxicircle.fa 		\
				--in_lib ${file}						\
				--out_prefix ${work_directory}/3-2nd_align/${file_name}.	
		done

		## 5. UMI check & merge
		if [ ${UMI} = T ]
		then
		python ${bin}/main/taf_UMI.py 	\
                        -i ${strain}.${gene}.${strand}	\
                        -w ${work_directory}
		fi
	done
done

if [ ${UMI} = F ]
then
for gene in $(ls ${Ref_gene})
do 
	for strand in fwd rev
	do 
		cat ${work_directory}/3-2nd_align/${strain}.${gene}.${strand}.*taf >> ${work_directory}/4-UMI_merge/taf/${strain}.${gene}.${strand}_dedup.taf
		cp ${work_directory}/2-fastq_extract/${strain}.${gene}.${strand}.fastq  ${work_directory}/4-UMI_merge/fastq/${strain}.${gene}.${strand}_dedup.fastq
	done
done
fi

## 6. select alignment
python ${bin}/main/select.py \
	-s ${strain} 	\
	-g ${Ref_gene}	\
	-w ${work_directory}	\
	-r ${RNA_reverse}







