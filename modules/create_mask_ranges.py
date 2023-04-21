from argparse import ArgumentParser, ArgumentTypeError
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import json
import pysam
from pathlib import Path


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')


def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument('--ref_dir',
                        type=str,
                        help='where the reference files are stored',
                        required=True)
    parser.add_argument('-o',
                        '--out_dirs',
                        nargs='+',
                        help='Output dir paths space separated to add multiple -d /path/one /path/two /path/three',
                        required=True)
    parser.add_argument('-d',
                        '--db_suffix',
                        help='suffix of the db i.e. "gen"',
                        required=False)
    parser.add_argument('--result_dir',
                        type=str,
                        help='where the config files are stored',
                        required=True)
    parser.add_argument('--prefix',
                        type=str,
                        default='',
                        help='optional identifying prefix if you need to save files to similar folder',
                        required=False)
    parser.add_argument('--threshold',
                        type=int,
                        default=27,
                        help='max depth of coverage threshold of set for exclusion mask',
                        required=False)




def get_fasta_seq_dict(fasta_path):
    fasta_seq_dict = {}
    genome_fasta_open = pysam.Fastafile(fasta_path)
    for ref in genome_fasta_open.references:
        fasta_seq_dict[ref] = genome_fasta_open.fetch(ref)
    genome_fasta_open.close()

    return fasta_seq_dict


def get_process_arrays_args():
    """	Inputs arguments from bash
    Gets the arguments, checks requirements, returns a dictionary of arguments
    Return: args - Arguments as a dictionary
    """
    parser = ArgumentParser()
    parse_process_arrays_args(parser)
    args = parser.parse_args()
    return args


args = get_process_arrays_args()
out_dirs = args.out_dirs
ref_dir = args.ref_dir
prefix = args.prefix
result_dir = args.result_dir
db_suffix = args.db_suffix
threshold = args.threshold
pd.options.mode.chained_assignment = None
os.makedirs(result_dir, exist_ok=True)
df = pd.DataFrame()

for out_dir in out_dirs:
    out_dir = os.path.join(out_dir, 'depth_all')
    file_list = os.listdir(out_dir)
    file_list = [os.path.join(out_dir, x) for x in file_list if
                 not x.startswith('._') and (f'.depth_all-{db_suffix}.csv' in x)]
    sample_list = []
    len_file_list = len(file_list)
    i = 0
    for file_i in file_list:
        i += 1
        sample_i = os.path.basename(file_i).split('.')[0]
        # chart progress
        print(f'{sample_i} - {i} of {len_file_list} - {out_dir}')
        # print(sample_i, i, len(baylor_id_list))
        # # sample_list.append(sample_i)
        # if not sample_i in baylor_id_list:
        #     continue
        # print(sample_i)
        df_i = pd.read_csv(file_i)
        df_i['SAMPLE'] = sample_i
        df_i['MCM_TYPE'] = [x.split('_')[0] for x in df_i['ALLELE']]
        df_i['ALLELE_CLASS'] = [x.split('_')[1] for x in df_i['ALLELE']]
        df_i['REF_LEN'] = df_i.groupby('ALLELE')['POSITION'].transform('max')
        # fix a typo in the header of a fasta (as needed)
        df_i['ALLELE_CLASS'] = [x if x != 'NKG2C-3-gen' else 'NKG2C3-gen' for x in df_i['ALLELE_CLASS']]
        df = pd.concat([df, df_i], ignore_index=True)

class_list = list(df['ALLELE_CLASS'].unique())
j = 0
len_class_list = len(class_list)
for class_i in class_list:
    j += 1

    print(f'{class_i} - {j} of {len_class_list}')
    # df[df['ALLELE_CLASS'] == class_i].to_csv(os.path.join(result_dir, f'{class_i}_{suffix}_depth.csv'), index=False)
    df[df['ALLELE_CLASS'] == class_i].to_csv(os.path.join(result_dir, f'{prefix}{class_i}_depth.csv.gz'), index=False)

# df_max = pd.DataFrame()
class_list = ['CD94-gen',
              'NKG2A-gen',
              'NKG2C1-gen',
              'NKG2C2-gen',
              'NKG2C3-gen',
              'NKG2D-gen',
              'NKG2F-gen']
# for class_i in class_list:
#     print(class_i)
#     df_i = pd.read_csv(os.path.join(result_dir, f'{class_i}_{suffix}_depth.csv.gz'))
#     df_max_i = df_i.groupby(['ALLELE','POSITION'])['DEPTH'].agg(max).reset_index()
#     df_max = pd.concat([df_max, df_max_i], ignore_index=True)
#     df_max_i.to_csv(os.path.join(result_dir, f'{class_i}_{suffix}_depth_max_allele.csv.gz'), index=False)
df_max = df.groupby(['ALLELE', 'POSITION'])['DEPTH'].agg(max).reset_index()
df_max.to_csv(os.path.join(result_dir, '{prefix}depth_max_allele.csv'), index=False)


# df_max = pd.read_csv(os.path.join(result_dir, '{suffix}_depth_max_allele.csv'))
df_max['REF_LEN'] = df_max.groupby('ALLELE')['POSITION'].transform('max')

allele_list = list(df_max['ALLELE'].unique())
df_range = pd.DataFrame()
# Array optimized tricks to determine the ranges
for allele_i in allele_list:
    df_max_i = df_max[df_max['ALLELE'] == allele_i]
    pos_max_1 = max(df_max_i['POSITION']) + 1
    df_max_i = pd.concat([pd.DataFrame({'ALLELE': [allele_i], 'POSITION': 0, 'DEPTH': threshold + 1}),
                          df_max_i,
                          pd.DataFrame({'ALLELE': [allele_i], 'POSITION': pos_max_1, 'DEPTH': threshold + 1})],
                         ignore_index=True)
    df_max_i = df_max_i[df_max_i['DEPTH'] >= threshold]
    pos_list = list(df_max_i['POSITION'])
    pos_list.pop(0)
    pos_list.append(max(pos_list) + 1)
    df_max_i['DIFF'] = pos_list - df_max_i['POSITION']
    # turn the diff into a range if it is greater than one
    df_range_i = df_max_i[df_max_i['DIFF'] > 1]
    df_range_i['START'] = df_range_i['POSITION'] + 1
    df_range_i['END'] = df_range_i['POSITION'] + df_range_i['DIFF'] - 1
    df_range_i = df_range_i[['ALLELE', 'START', 'END']]

    df_range = pd.concat([df_range, df_range_i], ignore_index=True)

df_range.to_csv(os.path.join(ref_dir, 'mask_range.csv'), index=False)
print(f'Outputted mask_range.csv to: {ref_dir}')
# Create a masked fasta to help with alignment bam files
# Import as a dictionary.
fasta_seq_dict = get_fasta_seq_dict(fasta_path=os.path.join(ref_dir, 'bait.fasta'))


def test_ranges(position, start_list, end_list):
    for start_i, end_i in zip(start_list, end_list):
        if (position >= start_i) and (position <= end_i):
            return True
    return False


with open(os.path.join(ref_dir, f'{prefix}nmask.fasta'), 'w') as f:
    # loop through the ranges and if the position is in the ranges (inclusive) replace the BP character with N
    for allele_i, seq_i in fasta_seq_dict.items():
        df_range_i = df_range[df_range['ALLELE'] == allele_i]
        start_list = list(df_range_i['START'])
        end_list = list(df_range_i['END'])
        if len(end_list) > 0:
            seq_n = ['N' if test_ranges(y, start_list, end_list) else x for x, y in
                     zip(seq_i, range(1, len(seq_i) + 1))]
            seq_str = ''.join(seq_n)
        else:
            seq_str = seq_i
        f.write(f'>{allele_i}\n')
        f.write(f'{seq_str}\n')
print(f'Outputted {prefix}nmask.fasta to: {ref_dir}')
