# Tless Alignment Toolkit
## Work environment
Python 3.9  
Linux CentOS7

### Packages
#### alignlib  
A tool within T-Aligner3.3 package used to align sequencing data to Tless reference. 
##### Github repository  
        https://github.com/jalgard/T-Aligner3.3 

#### python-Levenshtein  
A package used to calculate editing distance between two sequence.
##### Installation
        pip install python-Levenshtein
	
## Usage
    sh D_T-Aligner_main.sh <strain_name> <read_library_path> <ref_gene_path> <have_UMI (T/F)> <toolkits_bin> <output> <RNA_reverse (T/F)> 
    sh D_T-Aligner_summary.sh <strain_name> <path_to_taf> <ref_gene_path> <toolkits_bin> <output>
### 1. D_T-Aligner_main.sh

#### Input
##### (1) <strain_name>
The prefix of your fastq file.  
e.g. <strain name>.fastq
##### (2) <read_library_path>  
Pathway to your fastq file.  
The full path will be like:  
<read_library_path>/<strain_name>.fastq
##### (3) <ref_gene_path>  
Pathway to your reference gene that sequencing reads will be aligned to.  
*You may need to change the ref sequence’s path & filename to match the script:  
--in_ref <ref_gene_path>/<gene_name>/<gene_name>_maxicircle.fa  
Where <gene_name> is found by:  
for gene in $(ls <ref_gene_path>)
##### (4) <have_UMI (T/F)>  
Whether the sequencing library have UMI labeled in fastq file that will be used to remove duplication.  
T: Have; F: Don’t have.
##### (5) <toolkits_bin>
Directory to Tless alignment toolkits scripts bin.
##### (6)  output
Pathway for output.
##### (7) <RNA_reverse (T/F)>
Whether using your reverse complement sequencing read as input for RNA-seq.

#### Running toy eCLIP data against 20genes (12 cryptogenes and 8 never edit genes) in the playground attached (UMI=T, RNA_reverse=F):
	sh ./toolkits/D_T-Aligner_main.sh toy_eCLIP ./playground/fastq ./playground/20-genes T ./toolkits/bin ./playground/toy_eCLIP F    
#### Running toy CO3 data against co3 gene in the playground attached (UMI=T, RNA_reverse=F):
	sh ./toolkits/D_T-Aligner_main.sh toy_CO3 ./playground/fastq ./playground/co3_gene T ./toolkits/bin ./playground/toy_CO3 F    



#### Output
##### (1) 0-split_fastq  
Store the split fastq reads.
##### (2) 1-1st_align  
Store the alignlib output for the split fastq.
##### (3) 2-fastq_extract  
Store the stranded reads and the strandness-failed-to-determind reads as well as the split version.
##### (4) 3-2nd_align  
Store the alignlib output for the split stranded reads.
##### (5) 4-UMI_merge  
 _4-UMI_merge/fastq_  and  _4-UMI_merge/taf_  store the merged stranded read and UMI filter result (if have UMI).
##### (6) 5-select_final  
Final fastq and T-Aligner output after selection of one read that was aligned to multiple genes.  
Besides, the T-Aligner result will be converted into sam format. The CIGAR will be set as xxSyyMzzS, where yyM is the precise description. The remaining Softclip will be evenly assign to xxS and zzS. In the optional field, 'NM' refers to the sum of T that will be insert/delete in the read and 'EM' refers to editing matrix that produced by T-Aligner.




### 2. D_T-Aligner_summary.sh

#### Input
##### (1) <strain_name>
The prefix of your fastq file.  
##### (2) <path_to_taf>
Pathway to your taf file produced in  **1. D_T-Aligner_main.sh** .  
##### (3) <ref_gene_path>
Pathway to your reference gene that sequencing reads will be aligned to.   
##### (4) <toolkits_bin>
Directory to Tless alignment toolkits scripts bin.  
##### (5) output
Pathway for output.  

#### Running toy eCLIP data against 20genes (12 cryptogenes and 8 never edit genes) in the playground attached (UMI=T, RNA_reverse=F):
	sh ./toolkits/D_T-Aligner_summary.sh toy_eCLIP ./playground/toy_eCLIP/5-select_final/taf ./playground/20-genes ./toolkits/bin ./playground/toy_eCLIP/6-summary    
#### Running toy CO3 data against co3 gene in the playground attached (UMI=T, RNA_reverse=F):
	sh ./toolkits/D_T-Aligner_summary.sh toy_CO3 ./playground/toy_CO3/5-select_final/taf ./playground/co3_gene ./toolkits/bin ./playground/toy_CO3/6-summary  

#### Output  
##### (1) 6-summary  
There will be two kind of summary, read count and editing plot.  
Within  _6-summary/reads_count_ , the fwd/rev alignment count and the type of the alignment for each reference gene will be given. For each strandness, the order will be <gene_name> <total_count> <pre_edit_count> <partially_edited_count> <fully_edited_count>.  
Within  _6-summary/plot_ , editing plot in pdf, tiff, compressed tiff format will be given. For plot in pdf format, the plot label is enabled.

	

