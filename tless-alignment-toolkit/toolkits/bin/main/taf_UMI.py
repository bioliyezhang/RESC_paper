import argparse, string
import glob

args=argparse.ArgumentParser()
args.add_argument('-i',type=str,help='prefix')
args.add_argument('-w',type=str,help='work path')
args=args.parse_args()
prefix='%s'%args.i
work_path='%s'%args.w

def get_dict_key(dic, value):
    key = list(dic.keys())[list(dic.values()).index(value)]
    return key

def read_taf(filename):
    seq_index=[]
    line_list=[]
    info=[]
    with open(filename, 'r') as ifile:
        for line in ifile:
            seq_index.append(int(line.split()[1][1:]))
            line_list.append(line.strip())
            info.append('\t'.join(line.split('\t')[2:5]))
    table_index=zip(seq_index,line_list,info)
        
    sorted_table_index = sorted(table_index,key=lambda x:x[0])
    seq_index,sorted_linelist,sorted_info=[list(x) for x in zip(*sorted_table_index)]
    return sorted_linelist,sorted_info

def read_fastq(filename):
    table=[]
    UMI_list=[]
    with open(filename,'r') as ifile:
        for line in ifile.readlines():
            if line[0]=='@':
                table.append([line.strip()])
                UMI_list.append(line.strip().split('_')[1])
            else:
                table[-1].append(line.strip())
    return table,UMI_list

def fastq_select(seqs,subgroup):
    if len(seqs) ==1:
        return subgroup[0],True
    length_table = [len(ele) for ele in seqs]
    seq_letter_table = [[0 for ele in range(max(length_table))] for ele in range(4)]
    match_table={
        'A':0,
        'T':1,
        'C':2,
        'G':3,
    }
    for seq_id in range(len(seqs)):
        for letter_id in range(len(seqs[seq_id])):
            seq_letter_table[match_table[seqs[seq_id][letter_id]]][letter_id]+=1
    count_table=[list(j) for j in zip(*seq_letter_table)]
    final_seq=''
    for posi in count_table:
        maxi=max(posi)
        if maxi< len(seqs)/2:
            break
        if posi.count(maxi) !=1:
            return subgroup,False
        letter_index=posi.index(maxi)
        final_seq+=get_dict_key(match_table,letter_index)
    if final_seq in seqs:
        return subgroup[seqs.index(final_seq)],True
    else:
        return subgroup,False

def UMI_check(UMI_list,sorted_linelist,sorted_info,fastq_table):
    index_t=[i for i in range(len(UMI_list))]
    info_table=zip(UMI_list,sorted_info,sorted_linelist,index_t)
    sorted_info_table = sorted(info_table,key=lambda x:x[0])
    UMI_final=[]
    info_final=[]
    line_final=[]
    fqindex_final=[]
    for i in range(len(sorted_info_table)):
        ele = sorted_info_table[i]
        UMI,info,line,fq_index=ele
        if UMI not in UMI_final:
            UMI_final.append(UMI)
            info_final.append([info])
            line_final.append([[line]])
            fqindex_final.append([[fq_index]])
        else:
            if info not in info_final[-1]:
                info_final[-1].append(info)
                line_final[-1].append([line])
                fqindex_final[-1].append([fq_index])
            else:
                index = info_final[-1].index(info)
                line_final[-1][index].append(line)
                fqindex_final[-1][index].append(fq_index)
    # duplication check
    with open(work_path+'/4-UMI_merge/fastq/'+prefix+'_dedup.fastq','a') as fq_out:
        with open(work_path+'/4-UMI_merge/fastq/fail/'+prefix+'_fail_dedup.fastq','a') as fq_f_out:
            with open(work_path+'/4-UMI_merge/taf/'+prefix+'_dedup.taf','a') as taf_out:
                for i in range(len(UMI_final)):
                    for subgroup in fqindex_final[i]:                                   
                        seqs = [fastq_table[ele][1] for ele in subgroup]
                        index,state = fastq_select(seqs, subgroup)
                        if state:
                            for ele in fastq_table[index]:
                                fq_out.write(ele+'\n')
                            taf_out.write(sorted_linelist[index]+'\n')
                        else:
                            for id in index:
                                for ele in fastq_table[id]:
                                    fq_f_out.write(ele+'\n')
                            fq_f_out.write('\n')



if __name__ == '__main__':
    t_sorted_linelist,t_sorted_info,t_fastq_table,t_UMI_list=[],[],[],[]
    suffix_list=[]
    for filename in glob.glob(work_path+'/2-fastq_extract/split/'+prefix+'*'): 
        suffix_list.append(filename.split('fastq.')[1])

    for suffix in suffix_list:
        fastq = work_path+'/2-fastq_extract/split/'+prefix+'.fastq.'+suffix
        taf_in = work_path+'/3-2nd_align/'+prefix+'.fastq.'+suffix+'.mapped_reads.taf'
        sorted_linelist,sorted_info=read_taf(taf_in)
        fastq_table,UMI_list=read_fastq(fastq)
        
        t_sorted_linelist += sorted_linelist
        t_sorted_info += sorted_info
        t_fastq_table += fastq_table
        t_UMI_list += UMI_list
    UMI_check(t_UMI_list,t_sorted_linelist,t_sorted_info,t_fastq_table)


