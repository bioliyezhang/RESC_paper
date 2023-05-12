import  matplotlib,os,sys
from matplotlib import pyplot as plt
import argparse
matplotlib.use('Agg')

def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('-i', 
        action='store',
        help="taf path")

    optionParser.add_argument('-s', 
        action='store',
        help="strain name")

    optionParser.add_argument('-g',
        action='store',
        help="Ref Gene pathway")

    optionParser.add_argument('-o',
        action='store',
        help="out path")


    return optionParser

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

def ReadName(input_lines):
    name = []
    for line in input_lines:
        name.append(line.strip())
    return name

def multi_mkdir(path):
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

if __name__ == '__main__':
    runArgs = GetOptParser().parse_args(sys.argv[1:])
    
    font = {'family': 'Arial',                                                                      # output font
         'weight': 'normal',
         'size': 20,
         }

    strain=runArgs.s               
    gene_path=runArgs.g
    in_path=runArgs.i
    out_path=runArgs.o
    gene_list=[]
    for gene_name in os.listdir(gene_path):
        gene_list.append(gene_name)
    gene_list=sorted(list(set(gene_list)))                                                              
    plt.figure(figsize=(6.4*1.5,4.8*1.5)) 
    for gene in gene_list:
        for strand in ['fwd','rev']:
            inputTafFile = in_path+'/'+strain+'.'+gene+'.'+strand+'.final.taf'          
            inputFastaData = open(gene_path+'/'+gene+'/'+gene+'_maxicircle.fa', 'r').readlines()        
            fastaEntries = ReadFasta(inputFastaData)
            inputFastaName = 0
            if None is not None:
                for i, n in enumerate(fastaEntries):
                    if n[0] == None:
                        inputFastaName = i
                        break

            referenceSequence = fastaEntries[inputFastaName][1]
            # get cryptogene reference length
            ref_tless_length = 0

            for i in referenceSequence.lstrip().rstrip().upper():
                if i != 'T':
                    ref_tless_length += 1

            # fill matrices with zero values
            es_read_matrix = []
            es_productive_matrix = []

            for y in range(-20, 20):
                es_read_matrix.append([0 for x in range(ref_tless_length)])
                es_productive_matrix.append([-0.01 for x in range(ref_tless_length)])

            # read TAF

            


            with open(inputTafFile, 'r') as ifile:
                for line in ifile:
                    toks = line.split('\t')
                    s = int(toks[2])
                    es = [int(x) for x in toks[4][:-1].split(';')[1:-1]]
                    # skip VERY long insertions
                    skip = False
                    for ess in es:
                        if ess > 20:
                            skip = True
                            break
                    if skip:
                        continue

                    # detect if editing drops to reference (cutoff 7)
                    for p, ess in enumerate(es):
                        if s + p >= ref_tless_length:
                            continue
                        es_read_matrix[20-ess][s + p] += 1
                        if p > 7:
                            if es[p-8:p] == [0, 0, 0, 0, 0, 0, 0, 0]:
                                es_productive_matrix[20-ess][s + p] += 1

            l_from = 0
            l_to   = ref_tless_length
            if None is not None:
                f, t = str(None).split(',')
                l_from = int(f)
                l_to   = int(t)
            
            covi = [0.0 for x in range(len(es_read_matrix[0][l_from:l_to]))]
            inserti = [0.0 for x in range(len(covi))]
            deletioni = [0.0 for x in range(len(covi))]
            xi   = [x for x in range(len(covi))]
            
            for y in range(len(es_read_matrix)):
                for x in range(len(covi)):
                    if y == 20:
                        covi[x] += es_read_matrix[y][x]
                    elif y<20:
                        inserti[x] += es_read_matrix[y][x]
                    elif y>20:
                        deletioni[x] += es_read_matrix[y][x]
            str_inserti = [str(ele) for ele in inserti]
            str_deletion=[str(ele) for ele in deletioni]
            str_covi = [str(ele) for ele in covi]
            p1 = plt.bar(xi, covi,label='unedit',width=1.0)
            p2 = plt.bar(xi, inserti,bottom=covi,label='insertion',width=1.0)
            p3 = plt.bar(xi, deletioni,bottom=[covi[i]+inserti[i] for i in range(len(covi))],label='deletion',width=1.0)

            plt.xlabel('Tless position',font,y=-1.02)
            plt.ylabel('Reads count',font,x=-1.02)
            plt.xticks(fontproperties = 'Arial', size = 20)
            plt.yticks(fontproperties = 'Arial', size = 20)
            plt.subplots_adjust(left=0.15, bottom=0.15)
            ax=plt.gca()
            ax.spines['bottom'].set_linewidth(1.5)
            ax.spines['left'].set_linewidth(1.5)
            ax.spines['right'].set_linewidth(1.5)
            ax.spines['top'].set_linewidth(1.5)
            output_pic = out_path+'/tif/bars.'+strain+'.'+gene+'.'+strand+'.tif'
            plt.savefig(output_pic,dpi=300,bbox_inches='tight')
            plt.title(strain.split('_')[0]+'_'+gene+'_'+strand,font,y=1.02) 
            plt.legend()
            output_pic = out_path+'/pdf/bars.'+strain+'.'+gene+'.'+strand+'.pdf'
            plt.savefig(output_pic,dpi=300,bbox_inches='tight')
            plt.clf()

    