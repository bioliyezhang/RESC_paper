#python3
import argparse
import sys
import re
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import ttest_rel
from scipy.stats import ttest_ind
from scipy.stats import wilcoxon
from statannotations.Annotator import Annotator

parser = argparse.ArgumentParser()
parser.add_argument("-i1", dest="i1")
args = parser.parse_args()

dat = pd.read_csv(args.i1,sep='\t')

fig,ax = plt.subplots()
#pal={"RESC1":"#37BC9B","RESC2":"#aa5500","RESC5":"FC6E51","RESC6":"5D9CEC"}
pal={"RESC1":(0.216,0.737,0.608),"RESC2":(0.667,0.333,0),"RESC5":(0.988,0.431,0.318),"RESC6":(0.365,0.612,0.925)}

ax = sns.violinplot(x="sample",y="Utail_length",data=dat,palette=pal)

df1=dat[dat['sample']=="RESC1"]["Utail_length"]
df2=dat[dat['sample']=="RESC2"]["Utail_length"]
df5=dat[dat['sample']=="RESC5"]["Utail_length"]
df6=dat[dat['sample']=="RESC6"]["Utail_length"]

box_pairs = [("RESC2","RESC5"),("RESC5","RESC6"),("RESC1","RESC5"),("RESC1","RESC6"),("RESC2","RESC6")]
annotator =  Annotator(ax, data=dat, x="sample",y="Utail_length",pairs=box_pairs)

annotator.configure(test='Mann-Whitney', text_format='star', comparisons_correction="fdr_bh", line_height=0.03,  line_width=1)
annotator.apply_test().annotate(line_offset_to_group=0.2, line_offset=0.1)

plt.xlabel("")
plt.savefig(args.i1.split(".")[0]+'.'+'pdf')