import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='filter primers for each oligo group')
parser.add_argument('input_file', help='Path to input primer files')
parser.add_argument('output_file', help='Path to output file')
args = parser.parse_args()

'''
The input file is required to have this format
and this file can be generated from running 'get_all_primers.py' script
an example of this file is 19gbks_output3_test_all_primers.txt at
/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/test2/T3Pio/T3Pio_Main
OG0003142primerGroup0   ATGAAACCGTATCTGTCCCGC   TRGGACGAATCACCGGAGCTG
OG0003142primerGroup2   GCGGAACAGGACGTAAAAGAC   GTCTTCAGCGCCAAAGGTTTT
OG0003142primerGroup3   GCGGAACAGGACGTAAAAGAC   GTTTTCGCCAGTTCAATCGCT
OG0003142primerGroup8   GCGGAACAGGACGTAAAAGAC   TTCAATCGCTTCGTTGATGGC
'''
df = pd.read_csv(args.input_file, sep='\t', header=None, names=['Group', 'forward', 'reverse'])
df['group_number'] = df['Group'].str.split('primerGroup').str[0]

#Primer scoring block. Weights the score based on position of variable site(3' end bad) and length of primer so longer primers with a variable base in the same 5' position as a shorted primer will be be prefered. 
#Both left and right primer scores are combined for the total score for comparison
def primerScore(primer):
    spot = 0
    position = 1
    for letter in primer:
        if letter not in ('ATGC'):
            spot = position
        position += 1
    score = spot/(len(primer))
    return score

def get_lowest_score(group):
    scores = group.apply(lambda row: primerScore(row["forward"]) + primerScore(row["reverse"]), axis=1)
    lowest_score = scores.min()
    lowest_score_rows = group[scores == lowest_score]
    return lowest_score_rows.iloc[0]  # Select the first row with the lowest score if there is a tie

lowest_score_rows = df.groupby("group_number").apply(get_lowest_score).reset_index(drop=True)

lowest_score_rows = lowest_score_rows[['Group', 'forward', 'reverse']]
lowest_score_rows.to_csv(args.output_file, sep="\t", header=False, index=False)


