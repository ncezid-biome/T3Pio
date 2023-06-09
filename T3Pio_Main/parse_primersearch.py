import re
import argparse
import glob

'''
this script parse the primersearch result and creates a dictionary with
primer name as the key and a list of contig name as the value

this script is intended to be run in folder:
/scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/test2/T3Pio/T3Pio_Main

And can be run as:
python3 parse_primersearch.py /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/StoolBugsMultifastas_primersearch 
            /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/T3Pio_Data/StoolBugsMultifastas_primersearch/combined_kraken123_no_Salmonella_contigs

'''
parser = argparse.ArgumentParser(description='parse primersearch result')
parser.add_argument('input_dir', help='Path to input directory of .ps files (primersearch)')
# parser.add_argument('output_file', help='Path to output file')
parser.add_argument('bad_contigs', help='Path to bad contigs list file')
args = parser.parse_args()

def merge_dictionaries(dict1, dict2):
    merged_dict = dict1.copy()

    for key, value in dict2.items():
        if key in merged_dict:
            merged_dict[key].extend(value)
        else:
            merged_dict[key] = value

    return merged_dict


def parse_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = {}
    current_primer = None

    for line in lines:
        line = line.strip()
        if line.startswith('Primer name'):
            current_primer = line.split(' ')[-1]
            data[current_primer] = []
        elif line.startswith('Sequence:'):
            sequence = re.search(r'Sequence: (.+)', line).group(1)
            data[current_primer].append(sequence.strip())

    return data

parsed_data = None
for input_file in glob.glob(f"{args.input_dir}/*.ps"):
    if parsed_data is None:
        parsed_data = parse_file(input_file)
    else:
        parsed_data = merge_dictionaries(parsed_data, parse_file(input_file))

print (len(parsed_data))

#remove empty values
parsed_data = {primer: amplimers for primer, amplimers in parsed_data.items() if amplimers}

print (len(parsed_data))

with open (args.bad_contigs, 'r') as file:
    bad_contigs_list = [line.strip() for line in file.readlines()]

# Iterate over the dictionary and remove key/value pairs
parsed_data_filtered = {key: value for key, value in parsed_data.items() if not any(amplimer in bad_contigs_list for amplimer in value)}
## the equivalence of the above line if we use the 'good_contigs_list' instead of the 'bad_contigs_list'
# parsed_data_filtered = {key: value for key, value in parsed_data.items() if all(amplimer in bad_contigs_list for amplimer in value)}


print (len(parsed_data_filtered))

# Accessing the parsed data
# for primer, amplimers in parsed_data.items():
#     print(f"Primer: {primer}")
#     for amplimer in amplimers:
#         print(f"Amplimer: {amplimer}")
#     print()
    
    
# import random
# # Get 5 random keys from the dictionary
# random_keys = random.sample(parsed_data_filtered.keys(), 5)
# for key in random_keys:
#     print(f"Primer: {key}")
#     for amplimer in parsed_data[key]:
#         print(f"Amplimer: {amplimer}")
#     print()


with open ('19gbks_output3_test_all_primers_list', 'r') as file:
    primer_list = [line.strip() for line in file.readlines()]
    
print (f"original primer list is of size: {len(primer_list)}")

final_primer_list = list(set(parsed_data_filtered.keys()) & set(primer_list))

print (f"final primer list is of size: {len(final_primer_list)}")

#get only the oligo name part from the primer name
final_og_list = [primer_name.split('primer')[0] for primer_name in final_primer_list]

with open ('2461_primers_list', 'r') as file:
    primer_2461_list = [line.strip() for line in file.readlines()]
   
og_2461_list = [primer_name.split('primer')[0] for primer_name in primer_2461_list] 

print (f"total common oligo group count is: {len(set(final_og_list) & set(og_2461_list))}")

final_final_primer_list = list(set(final_primer_list) & set(primer_2461_list))
print (f"total common primer list is of size: {len(final_final_primer_list)}")
print (f"total common primer list before this filtering is: {len(set(primer_2461_list) & set(parsed_data.keys()))}")

og_3314_list = [primer_name.split('primer')[0] for primer_name in parsed_data.keys()]
print (f"total common oligo group count before this filtering is: {len(set(og_3314_list) & set(og_2461_list))}")
print (f"the missing og is: {(set(og_2461_list) - set(og_3314_list))}")
    


# extra_primer = list(set(primer_2461_list) - set(final_primer_list))
# print (f"extra primer size is: {len(extra_primer)}")
# random_keys = random.sample(extra_primer, 5)
# for key in random_keys:
#     print(f"extra Primer ?: {key}")