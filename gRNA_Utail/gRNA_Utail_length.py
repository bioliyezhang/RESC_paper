#python3
import sys
import matplotlib
import matplotlib.pyplot as plt
import pysam
matplotlib.use('Agg')
def parse_CIGAR_Softclip(CIGAR):
    start=''
    for ele in CIGAR:
        if ele.isdigit():
            start+=ele
        else:
            if ele=='S':
                return 'SoftClip',int(start)
            else:
                return 'Not',0

def count_element(list):
    result={}
    for i in set(list):
        result[i]=list.count(i)
    return result

bam = pysam.AlignmentFile(sys.argv[1],'rb')

soft_clip_size_list=[]
for i in bam:
    CIGAR=i.cigar
    seq=i.query_sequence    
    if i.cigar[-1][0]==4:
        ex_seq=seq[-CIGAR[-1][1]:]        
        if ex_seq.count('T')/len(ex_seq)==1.0:
            soft_clip_size_list.append(i.cigar[-1][1])
if len(soft_clip_size_list)==0:    
    print(str(0))
else:
    res=count_element(soft_clip_size_list)
    x=[]
    y=[]
    for k,v in res.items():
        x.append(k*v)
        y.append(v)
    total1=sum(x)
    total2=sum(y)
    ave=total1/total2
    print(str(ave))