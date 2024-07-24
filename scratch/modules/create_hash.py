from argparse import ArgumentParser, ArgumentTypeError
import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import math
import random


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
    parser.add_argument('--bait_fasta',
                        type=str,
                        help='Where you ipd-miseq database combined fasta exists',
                        default=None,
                        required=False)
    parser.add_argument('--ref_dir',
                        type=str,
                        help='where the config files are stored',
                        required=True)
    parser.add_argument('--read_length',
                        type=int,
                        help='Minimum read length i.e. 151. If the reads are trimmed you may need shorter 151.',
                        default=151,
                        required=False)
    parser.add_argument('--file_split',
                        type=int,
                        help='how many df to split hash into to reduce ram',
                        default=10,
                        required=False)


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

bait_fasta = args.bait_fasta
ref_dir = args.ref_dir
read_length = args.read_length
file_split = args.file_split
pd.options.mode.chained_assignment = None

if (ref_dir is None) and (bait_fasta is None):
    print("A ref_dir or  bait_fast  must be declared")
    exit()
if bait_fasta is None:
    bait_fasta = os.path.join(ref_dir, 'bait.fasta')


def fasta_to_df(fasta_path=None, header_name='ALLELE', sequence_name='SEQUENCE', as_df=False):
    fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
    fasta_dict = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_dict[name] = sequence
    fasta_sequences.close()
    if as_df:
        return pd.DataFrame(fasta_dict.items(), columns=[header_name, sequence_name])
    return fasta_dict


def make_lists(seq, read_length=151):
    seq_l = len(seq)
    seq_list = []
    prev_seq = ''
    read_half = math.floor(read_length / 2)
    for x in range(0, seq_l - read_half):
        y = x + read_length
        # if seq[x:y] != prev_seq:
        seq_list.append(seq[x:y])
        # prev_seq = seq[x:y]
    return seq_list


def compare_full(df_exp_c, df_e,  file_count=10):
    # shoten the columns to save space and change type from float to int
    df_exp_compare_full = df_exp_c[['ALLELE_2', 'SEQ_LIST_2', 'START_2', 'END_2']]
    df_exp_compare_full.rename(columns={'SEQ_LIST_2': 'SEQ_LIST'}, inplace=True)
    print(df_exp_compare_full)
    print(df_exp_compare_full.columns)
    df_exp_compare_full['START_2'] = df_exp_compare_full['START_2'].astype(int)
    df_exp_compare_full['END_2'] = df_exp_compare_full['END_2'].astype(int)
    df_exp_full = df_e[['ALLELE', 'SEQ_LIST', 'START', 'END']]
    df_exp_full['START'] = df_exp_full['START'].astype(int)
    df_exp_full['END'] = df_exp_full['END'].astype(int)
    # df_exp_compare_full
    unique_allele_list = list(df_exp_compare_full['ALLELE_2'].unique())
    random.seed(4)
    random.shuffle(unique_allele_list)
    df_all = pd.DataFrame()
    for idx, chunk in enumerate(np.array_split(unique_allele_list, file_count)):
        df_chunk = df_exp_compare_full[df_exp_compare_full['ALLELE_2'].isin(chunk)]
        print(idx, file_count, len(df_chunk))
        df_i = df_exp_full.merge(df_chunk, on=['SEQ_LIST'], how='inner')
        # df_exp_full_merge = df_exp_full_merge[df_exp_full_merge['allele'] != df_exp_full_merge['allele_2']]
        df_i.drop(columns=['SEQ_LIST'], inplace=True)
        df_i.drop_duplicates(inplace=True)

        print('Sort values')
        df_i.sort_values(by=['ALLELE', 'ALLELE_2', 'START', 'END'], inplace=True)
        allele_2_list = list(df_i['ALLELE_2'])
        allele_2_list.insert(0, 'XZAS')
        allele_2_list.pop(-1)
        allele_list = list(df_i['ALLELE'])
        allele_list.insert(0, 'XZAS')
        allele_list.pop(-1)
        df_i['NEW_GROUP'] = [1 if ((x != y) or (z != w)) else 0 for x, y, z, w in zip(df_i['ALLELE_2'],
                                                                                      allele_2_list,
                                                                                      df_i['ALLELE'],
                                                                                      allele_list)]
        for col_i in ['START', 'START_2', 'END', 'END_2']:
            print(col_i)
            start_list = list(df_i[col_i])
            start_list.insert(0, 0)
            start_list.pop(-1)
            df_i['{0}_DIFF'.format(col_i)] = df_i[col_i] - start_list
            df_i['{0}_DIFF'.format(col_i)] = [np.nan if x > 0 else y for x, y in zip(df_i['NEW_GROUP'],
                                                                                     df_i['{0}_DIFF'.format(col_i)])]

        df_i['SUM_POS'] = df_i['START'] + df_i['END']
        print('Find max')
        df_i['MAX_POS'] = df_i.groupby(['ALLELE', 'ALLELE_2'])['SUM_POS'].transform('max')
        print(df_i)
        df_m = df_i[df_i['START_DIFF'].isna() | (df_i['START_DIFF'] > 1) | (df_i['END_DIFF'] > 1) | (
                    df_i['MAX_POS'] == df_i['SUM_POS'])]
        df_m['KEEP_LAST'] = [1 if ((x == y) and (not (z < 2))) else 0 for x, y, z in
                             zip(df_m['SUM_POS'], df_m['MAX_POS'], df_m['END_DIFF'])]
        df_m['END_SHIFT'] = [x if y < 2 else x - y for x, y in zip(df_m['END'], df_m['END_DIFF'])]
        df_m['LAST'] = [1 if (x == y) else 0 for x, y in zip(df_m['SUM_POS'], df_m['MAX_POS'])]
        df_m['END_2_SHIFT'] = [x if y < 2 else x - z for x, y, z in
                               zip(df_m['END_2'], df_m['END_DIFF'], df_m['END_2_DIFF'])]

        end_shift = list(df_m['END_SHIFT'])
        end_shift.pop(0)
        end_shift.append(np.nan)
        df_m['END_SHIFTED'] = end_shift

        end_shift = list(df_m['END_2_SHIFT'])
        end_shift.pop(0)
        end_shift.append(np.nan)

        df_m['END_2_SHIFTED'] = end_shift
        df_n = df_m[(df_m['KEEP_LAST'] > 0) | (df_m['LAST'] < 1)]
        df_n['END_2_SHIFTED'] = df_n['END_2_SHIFTED'].where(pd.notnull, df_n['END_2'])
        df_n['END_SHIFTED'] = df_n['END_SHIFTED'].where(pd.notnull, df_n['END'])
        df_n = df_n[['ALLELE', 'START', 'END_SHIFTED', 'ALLELE_2', 'START_2', 'END_2_SHIFTED']]
        df_n.rename(columns={'END_SHIFTED': 'END', 'END_2_SHIFTED': 'END_2'}, inplace=True)
        df_all = pd.concat([df_all, df_n], ignore_index=True)

    df_all.to_csv(os.path.join(ref_dir, 'ipd_hash.csv.gz'), index=False)

    del chunk
    del df_exp_compare_full
    del df_exp_full
    return

# Read the fasta as a dataframe
df = fasta_to_df(fasta_path=bait_fasta, header_name='ALLELE', sequence_name='SEQUENCE', as_df=True)
df['REF_SEQ_LEN'] = df['SEQUENCE'].str.len()
df['SEQ_LIST'] = df['SEQUENCE'].apply(make_lists, read_length=read_length)
df.drop(columns=['SEQUENCE'], inplace=True)
# create a sequence of every possibile segment of the reference sequence from 76 to 151 in length
df_exp = df.explode('SEQ_LIST')
df_exp['LEN'] = df_exp['SEQ_LIST'].str.len()
df_exp = df_exp.reset_index().rename(columns={'index': 'allele_num'})
df_exp['rank'] = df_exp.groupby('ALLELE')['allele_num'].rank('first')
df_exp['rank'] = df_exp['rank'].astype(int)

# needs to consider ref length
# first do an else statement
df_exp["START"] = df_exp["rank"] - 1
# half_len = math.ceil(read_length/2)

df_exp['END'] = df_exp['START'] + df_exp['LEN']

df_exp = df_exp[['ALLELE',
                 'REF_SEQ_LEN',
                 'SEQ_LIST',
                 'START',
                 'END',
                 'LEN']].drop_duplicates(keep='first')

df_exp = df_exp.reset_index()
print('create a copy of the fasta converted and expanded dataframe to compare')
df_exp_compare = df_exp.copy()
# make a dataframe that we will join to later for comparisons

df_exp_compare.rename(columns={'index': 'index_2',
                               'ALLELE': 'ALLELE_2',
                               'REF_SEQ_LEN': 'REF_SEQ_LEN_2',
                               'SEQ_LIST': 'SEQ_LIST_2',
                               'START': 'START_2',
                               'END': 'END_2',
                               'LEN': 'LEN_2'}, inplace=True)

df_exp.rename(columns={'index': 'index_1'}, inplace=True)

print('create compare full lengths and export to files .gz')
print(df_exp_compare)
compare_full(df_exp_compare, df_exp, file_count=file_split)
print('collecting excess variables')
# n = gc.collect()
print("Number of unreachable objects collected by GC:")
