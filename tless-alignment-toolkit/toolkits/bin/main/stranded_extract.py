import argparse, string
import Levenshtein

args=argparse.ArgumentParser()
args.add_argument('-i',type=str,help='input')
args.add_argument('-r',type=str,help='fastq')
args.add_argument('-o',type=str,help='out_prefix')
args.add_argument('-R',type=str,help='RNA_reverse',default='F')
args=args.parse_args()
taf_in='%s'%args.i
fastq='%s'%args.r
out_path='%s'%args.o
RNA_reverse='%s'%args.R


def DNA_Compliment(sequence):
    transtab=str.maketrans('ATCGatcg','TAGCtagc')
    string=sequence.translate(transtab)
    return string.upper()[::-1]
            
def Read_taf(filename):
    seq_index=[]
    filter_index=[]
    seq_table=[]
    with open(filename, 'r') as ifile:
        for line in ifile:
            if GC_check(line.split()[5]):
                seq_index.append((int(line.split()[1][1:])-1)*4+1)
                seq_table.append(line.split()[-3])
            else:
                seq_index.append((int(line.split()[1][1:])-1)*4+1)
                filter_index.append((int(line.split()[1][1:])-1)*4+1)
                seq_table.append(line.split()[-3])

    return seq_index,filter_index,seq_table


def GC_check(seq):
    if (seq.count('G')+seq.count('C'))/(seq.count('G')+seq.count('C')+seq.count('A'))<0.2:
        return False
    else:
        return True




if __name__ == '__main__':    
    f_filter=open(out_path+'.filtered.fastq','a')
    f_fwd=open(out_path+'.fwd.fastq','a')
    f_rev=open(out_path+'.rev.fastq','a')
    f_fail=open(out_path+'.fail.fastq','a')        
    
    seq_index,filter_index,seq_table =Read_taf(taf_in)          


    with open(fastq,'r') as f:
        index=0
        seq=[]
        stat=''
        for line in f.readlines():
            
            seq.append(line)
            if index % 4 == 1 and index in seq_index:
                stat='on'
                
            elif len(seq)!=0 and index%4 == 3 and stat == 'on':
                read_index=seq_index.index(index-2)
                seq_taf = seq_table[read_index]
                trim_seq_taf = seq_taf[8:-8]
                if (index-2) in filter_index:
                    for ele in seq:
                        f_filter.write(ele)
                elif Levenshtein.distance(seq_taf,seq[1])< Levenshtein.distance(seq_taf ,DNA_Compliment(seq[1])):
                    for ele in seq:
                        if RNA_reverse=='F':
                            f_fwd.write(ele)
                        else:
                            f_rev.write(ele)
                elif Levenshtein.distance(seq_taf,seq[1])> Levenshtein.distance(seq_taf ,DNA_Compliment(seq[1])):
                #elif string_similar(seq_taf,seq[1]) < string_similar(seq_taf ,DNA_Compliment(seq[1])):
                    for ele in seq:
                        if RNA_reverse=='F':
                            f_rev.write(ele)
                        else:
                            f_fwd.write(ele)
                else:
                    for ele in seq:
                        f_fail.write(ele)
                stat=''
            if len(seq)==4:     
                seq=[]              
            index+=1

