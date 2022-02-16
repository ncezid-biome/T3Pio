def blastnRun(genome,db,outfile):

    return_code = subprocess.check_output(['blastn','-db',db,'-query',genome,'-out',outfile,'-outfmt','6','-max_hsps','1','-best_hit_overhang','0.1','-best_hit_score_edge','0.1','-perc_identity','95'])


import argparse
import subprocess

parser = argparse.ArgumentParser(description='Runs blastn on a genome against your chosen database')

parser.add_argument('genome',type=str,help='File path to genome')

parser.add_argument('db',type=str,help='File path to blast db')

parser.add_argument('outFile',type=str,help='Name of the outfile')

args = parser.parse_args()

if __name__ == '__main__':

    blastnRun(args.genome,args.db,args.outFile)
