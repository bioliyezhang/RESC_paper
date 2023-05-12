import argparse, sys
from fileinput import filename


def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('-t',
        action='store',
        help="taf file")

    optionParser.add_argument('-m',
        action='store',
        help="editing matrix")

    return optionParser

def read_matrix(Filename):
    f_in = open(Filename,'r')
    matrix=[]
    for line in f_in.readlines():
        matrix.append([line.split()[0],line.split()[1].split(';')])
    return matrix

def Parse_taf(taf_path,editing_matrix):
    taf_f = open(taf_path,'r')
    pre_edit=0
    fully_edit=0
    partially_edit=0

    for line in taf_f.readlines():
        start,end,tok=line.split()[2:5]
        if tok.count('0')==len(tok)/2:
            pre_edit+=1
        else:
            taf_edit=tok[:-1].split(';')
            if editing_matrix[int(start):int(end)-1] == taf_edit:
                fully_edit+=1
            else:
                partially_edit+=1
    print(str(pre_edit+fully_edit+partially_edit)+'\t'+str(pre_edit)+'\t'+str(fully_edit)+'\t'+str(partially_edit))

if __name__ == '__main__':
    runArgs = GetOptParser().parse_args(sys.argv[1:])
    taf_path= runArgs.t
    matrix_file=runArgs.m
    file_name=taf_path.split('/')[-1]
    editing_matrix=read_matrix(matrix_file)
    ref_matrix=''
    for ele in editing_matrix:
        if ele[0] in file_name:
            ref_matrix=ele[1]
    if len(ref_matrix)==0:
        f_in=open(taf_path,'r')
        print(str(len(f_in.readlines()))+'\tNA\tNA\tNA')
    else:
        result=Parse_taf(taf_path,ref_matrix)

