import os
import argparse

# set up command-line argument parsing
parser = argparse.ArgumentParser()
parser.add_argument('dir_path', help='Path to directory containing primer files')
parser.add_argument('output_file', help='Path to output file')
args = parser.parse_args()

# set the directory path and output file path
dir_path = args.dir_path
output_file = args.output_file

# open the output file for writing
with open(output_file, 'w') as out_f:

    # loop through each file in the directory
    for filename in os.listdir(dir_path):
        if filename.endswith('Primers'):
            primer_group = filename.split('Primers')[0] + 'primerGroup'
            file_path = os.path.join(dir_path, filename)
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip().split('\t')
                    primer_num = line[0]
                    primers = '\t'.join(line[1:])
                    output_line = primer_group + primer_num + '\t' + primers
                    out_f.write(output_line + '\n')

