import rasterio as rio
import rasterio
import os,sys
import argparse
from tqdm import tqdm

def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('-i',
        action='store',
        help="in path")

    optionParser.add_argument('-s',
        action='store',
        help="strain name")

    optionParser.add_argument('-o',
        action='store',
        help="out path")

    return optionParser

if __name__ == '__main__':
    runArgs = GetOptParser().parse_args(sys.argv[1:])
    Input_path= runArgs.i

    Output_path = runArgs.o
    strain = runArgs.s
    pathDir= os.listdir(Input_path)

    for i in tqdm(range(len(pathDir))):
        
        rasterfile = Input_path+"/"+pathDir[i]
        if strain in rasterfile:
            rasterdata = rio.open(rasterfile)
            rasterdata2= rasterdata.read()

            profile = rasterdata.profile
            
            # compress method
            profile.update(
                compress='lzw',  #method: rle, lzw etc
                )
            out_put_name=Output_path +"/"+pathDir[i]

            with rasterio.open(out_put_name, mode='w', **profile) as dst:
                dst.write(rasterdata2)