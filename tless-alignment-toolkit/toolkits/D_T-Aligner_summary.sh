#!/bin/bash!

strain=$1
In_path=$2
Ref_gene=$3
bin=$4
work_path=$5

## 6. summary
mkdir -p ${work_path}/reads_count/tmp_${strain}
mkdir -p ${work_path}/plot/Editing/pdf
mkdir -p ${work_path}/plot/Editing/tif
# read count
RC_path=${work_path}/reads_count
for s in fwd rev
do	
	for i in $(ls ${In_path}/${strain}*${s}*taf)
	do
		python ${bin}/summary/check_taf.py \
			-t $i \
			-m ${bin}/summary/12genes_editing_matrix.txt \
			>> ${RC_path}/tmp_${strain}/${strain}.${s}.tmp
	done
done
for gene in $(ls ${Ref_gene})
do
	echo $gene >> ${work_path}/reads_count/tmp_${strain}/genelist
done
paste ${RC_path}/tmp_${strain}/genelist ${RC_path}/tmp_${strain}/${strain}.fwd.tmp ${RC_path}/tmp_${strain}/${strain}.rev.tmp > ${RC_path}/${strain}.readcounts.txt
rm -rf ${RC_path}/tmp_${strain}

# plot
python ${bin}/summary/taf_plot.py	\
	-i ${In_path}		\
	-s ${strain}		\
	-g ${Ref_gene}		\
	-o ${work_path}/plot/Editing
# compress
mkdir -p ${work_path}/plot/Editing/tif_compressed
python ${bin}/summary/compress.py	\
	-i ${work_path}/plot/Editing/tif	\
	-o ${work_path}/plot/Editing/tif_compressed	\
	-s ${strain}
