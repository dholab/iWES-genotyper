from Bio import SeqIO
import pandas as pd
import numpy as np
import os

from argparse import ArgumentParser, ArgumentTypeError


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
                        help='dir path to the bait.fasta that was created by thte workflow.',
                        required=True)
    parser.add_argument('--miseq_fasta',
                        type=str,
                        help='path to original miseq database file.',
                        required=True)
    parser.add_argument('--out_file',
                        type=str,
                        help='File path of the csv lookup dataframe output',
                        required=True)
    parser.add_argument('--db_type',
                        type=str,
                        default='miseq',
                        help='sample name that is prefix for file name of expanded maps such as: 104294',
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
miseq_fasta = args.miseq_fasta
bait_fasta = args.bait_fasta
out_file = args.out_file
db_type = args.db_type

# miseq_fasta = '/Users/dabaker3/github/iWES-genotyper/ref/Mamu/Mamu_MHC-miseq.fasta'
# bait_fasta = '/Volumes/T8/iwes/Mamu_class_1_2/bait.fasta'
# out_file = '/Volumes/T8/iwes/Mamu_class_1_2/bait_to_miseq2.csv'
# db_type = 'miseq'


def fasta_to_df(fasta_path=None, header_name='allele', sequence_name='SEQUENCE', as_df=False):
    fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
    fasta_dict = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_dict[name] = sequence
    fasta_sequences.close()
    if as_df:
        return pd.DataFrame(fasta_dict.items(), columns=[header_name, sequence_name])
    return fasta_dict


def like_join_haplo2(x, df_ref, column_i='SEQUENCE', name_i='allele'):
    name_list = []
    for idx, row in df_ref.iterrows():
        #if x in row[column_i]:
        if (row[column_i] in x) or (row[column_i] in rev_comp_st(x)):
            name_list.append(row[name_i])
    return name_list


def rev_comp_st(seq):
    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    seq = seq[::-1]
    return seq


df_miseq_ref = fasta_to_df(miseq_fasta, header_name='allele', sequence_name='SEQUENCE', as_df=True)
df_bait_ref = fasta_to_df(bait_fasta, header_name='allele', sequence_name='SEQUENCE', as_df=True, )
print(len(df_bait_ref))
df_chunk_array = np.array_split(df_bait_ref, 10)
df_bait = pd.DataFrame()
counter = 0
for df_i in df_chunk_array:
    print(counter)
    df_i[f'{db_type}_allele'] = [like_join_haplo2(x, df_miseq_ref, 'SEQUENCE', 'allele') for x in df_i['SEQUENCE']]
    df_bait = pd.concat([df_bait, df_i], ignore_index=True)
    counter += 1
df_bait.to_csv(out_file, index=False)
