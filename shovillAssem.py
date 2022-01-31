def shovillAssembler(outdir,R1,R2):
    
    return_code = subprocess.check_output(['shovill','--assembler','skesa','--trim','ON','--cpus','32','--ram','50.00','--outdir',outdir,'--R1',R1,'--R2',R2])
    


import argparse
import glob
import subprocess

parser = argparse.ArgumentParser(description='Runs Shovill on raw reads')

parser.add_argument('ReadFile',type=str,help='File containing base name of reads to be assembled')

parser.add_argument('OutDir',type=str,help='Path to the out directory')

args = parser.parse_args()


if __name__ == '__main__':

    with open(args.ReadFile,'r') as f:
        readList = f.readlines()
    f.close()

    for read in readList:
        read = read.rstrip()

        try:

            readpair =  sorted(glob.glob(read+'*'))

            readpairName = readpair[0].split('_')[0]
        
            outfile = args.OutDir+'/'+readpairName

            shovillAssembler(outfile,readpair[0],readpair[1])
        
        except IndexError:

            continue
