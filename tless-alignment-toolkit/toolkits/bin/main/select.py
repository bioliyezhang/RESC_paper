import argparse
import glob,Levenshtein



args=argparse.ArgumentParser(description="select best alignment")
args.add_argument('-s',type=str,help='strain')
args.add_argument('-g',type=str,help='geneRef path')
args.add_argument('-w',type=str,help='work path')
args.add_argument('-r',type=str,help='RNA_reverse')
args=args.parse_args()
rna_reverse='%s'%args.r
strain='%s'%args.s
path='%s'%args.g
work_path='%s'%args.w

def READ_all(strain,gene_list):
    fastq_all=[]
    taf_all=[]
    for gene in gene_list:
        for strand in ['fwd','rev']:
            with open(work_path+'/4-UMI_merge/taf/'+strain+'.'+gene+'.'+strand+'_dedup.taf','r') as taf_in:
                for line in taf_in.readlines():
                    taf_all.append([line.strip(),gene,strand])
            taf_in.close()
            with open(work_path+'/4-UMI_merge/fastq/'+strain+'.'+gene+'.'+strand+'_dedup.fastq','r') as fastq_in:
                count=0
                fastq_all.append([])
                for line in fastq_in.readlines():
                    if count == 4:
                        fastq_all.append([])
                        count=0
                    fastq_all[-1].append(line.strip())
                    count+=1
                
                
    all_table=sorted(list(zip(fastq_all,taf_all)),key=lambda x:x[0][0])
    return all_table        

def READ_genes(path,gene_list):
    gene_seq=[]
    for gene in gene_list:
        with open(path+'/'+gene+'/'+gene+'_maxicircle.fa','r') as f:
            for line in f.readlines():
                if line[0]=='>':
                    gene_seq.append([line.strip()[1:],''])
                else:
                    gene_seq[-1][-1]+=line.strip()
    return gene_seq

def DNA_Compliment(sequence):
    transtab=str.maketrans('ATCGatcg','TAGCtagc')
    string=sequence.translate(transtab)
    return string.upper()[::-1]

def select_alignment(table,ref_gene,gene_list,rna_reverse):
    seqs=[]
    for seq in table:
        if len(seqs)==0:
            seqs.append([seq[0][0],[seq]])
        elif seq[0][0] not in seqs[-1]:
            seqs.append([seq[0][0],[seq]])
        else:
            seqs[-1][-1].append(seq)
    final_seqs=[]
    for seq in seqs:
        if len(seq[1])==1:                                          # no multiple hits
            final_seqs.append([seq[0],seq[1][0]])
        else:
            match_length=[]
            for taf in seq[1]:
                start,end=taf[1][0].split('\t')[2:4]
                match_length.append(int(end)-int(start))
            max_id = match_length.index(max(match_length))
            if match_length.count(match_length[max_id]) == 1:       # Pick the longest Tless match
                final_seqs.append([seq[0],seq[1][max_id]])
            else:                                                   # Check Tless mismatch
                sub_list=[]
                for index in range(len(match_length)):
                    if match_length[index]== max(match_length):
                        sub_list.append(index)
                distance_list=[]
                fastq_Tless=seq[1][0][0][1].replace('T','')
                
                for ele in sub_list:
                    gene=seq[1][ele][1][1]
                    start,end=seq[1][ele][1][0].split('\t')[6:]
                    ref_seq=ref_gene[gene_list.index(gene)][1][int(start):int(end)]
                    if seq[1][ele][1][2]=='fwd' or rna_reverse=='T':
                        distance_list.append(Levenshtein.distance(fastq_Tless,ref_seq.replace('T','')))
                    else:
                        distance_list.append(Levenshtein.distance(fastq_Tless,DNA_Compliment(ref_seq).replace('T','')))
                min_list=[]
                for id in range(len(distance_list)):
                    if distance_list[id]== min(distance_list):
                        min_list.append(id)
                for id in min_list:
                    final_seqs.append([seq[0],seq[1][id]])      
    return final_seqs
    
def Write(final_seqs,strain,ref_gene,gene_list):
    for index in range(len(gene_list)):
        gene=gene_list[index]
        ref_name=ref_gene[index][0]
        ref_length=len(ref_gene[index][1])
        for strand in ['fwd','rev']:
            sam_out = open(work_path+'/5-select_final/sam/'+strain+'.'+gene+'.'+strand+'.final.sam','a')
            sam_out.write('@SQ\tSN:'+ref_name+'\tLN:'+str(ref_length)+'\n')
    for seq in final_seqs:
        gene = seq[1][1][1]
        strand=seq[1][1][2]
        taf_out=open(work_path+'/5-select_final/taf/'+strain+'.'+gene+'.'+strand+'.final.taf','a')
        fastq_out = open(work_path+'/5-select_final/fastq/'+strain+'.'+gene+'.'+strand+'.final.fastq','a')
        sam_out = open(work_path+'/5-select_final/sam/'+strain+'.'+gene+'.'+strand+'.final.sam','a')
        taf_out.write(seq[1][1][0]+'\n')
        for line in seq[1][0]:
            fastq_out.write(line+'\n')
        sam_out.write('\t'.join(write_sam(seq[1][1][0],seq[1][0],strand))+'\n')



def ReadFasta(input_lines, split_space = False):

    fasta_data = []

    for line in input_lines:
        line = line.rstrip()
        if len(line) < 1:
            continue
        if line[0] == '>':
            if split_space == False:
                fasta_data.append([line[1:],''])
            else:
                fasta_data.append([line[1:].split(' ')[0],''])
        else:
            fasta_data[-1][1] += line.lstrip().rstrip().upper()

    for i, entry in enumerate(fasta_data):
        fasta_data[i].append(len(entry[1]))

    return fasta_data

def ReadName(input_lines,suffix):
    name = []
    for line in input_lines:
        name.append(line.strip())
    return name

def write_sam(taf,fastq,strand):            # taf: one line, fastq: 4 lines matrix
    split_taf=taf.split()
    QNAME=fastq[0][1:]
    if strand=='fwd':
        FLAG='0'
        SEQ=fastq[1]
        QUAL=fastq[3]
    else:
        FLAG='16'
        SEQ=DNA_Compliment(fastq[1])
        QUAL=fastq[3][::-1]
    RNAME=split_taf[0]
    POS=str(int(split_taf[-2])+1)
    MAPQ='30'
    editing_distance=0
    indel=0
    for i in split_taf[4][:-1].split(';'):
        editing_distance+=abs(int(i))
        indel+=(int(i))
    Matched=len(split_taf[-3])
    SoftClip=len(SEQ)-Matched
    if SoftClip!=0:
        CIGAR=str(SoftClip//2)+'S'+str(Matched)+'M'+str(SoftClip-SoftClip//2)+'S'    
    else:
        CIGAR=str(Matched)+'M'
    RNEXT='*'
    PNEXT='0'
    TLEN='0'
    Optional_Fields='NM:i:'+str(editing_distance)+' EM:Z:'+split_taf[4]

    return [QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,RNEXT,PNEXT,TLEN,SEQ,QUAL,Optional_Fields]


    


if __name__ == '__main__':
    gene_list=[]
    for file_name in glob.glob(work_path+'/1-1st_align/'+strain+'*aa*taf'):
        gene_list.append(file_name.split('aa.')[1].split('.1st')[0])
    gene_list=sorted(list(set(gene_list)))
    ref_gene=READ_genes(path,gene_list)
    all_table=READ_all(strain,gene_list)
    final_seqs=select_alignment(all_table,ref_gene,gene_list,rna_reverse)
    Write(final_seqs,strain,ref_gene,gene_list)
