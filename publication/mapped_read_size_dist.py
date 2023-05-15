import pandas as pd
import numpy as np
import math
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
    parser.add_argument('--in_dir',
                        type=str,
                        help='filepath to the fasta. The fasta seq_ids must match ids in the  protein groups path',
                        required=True)
    parser.add_argument('--out_dir',
                        type=str,
                        help='path to ypur depth files.',
                        required=True)



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
# same arguments to a local variable by same name as the argument

out_dir = args.out_dir
in_dir = args.in_dir

filelist = os.listdir(in_dir)
filelist = [os.path.join(in_dir, x) for x in filelist if
            not x.startswith('._') and x.endswith('_expanded_maps-gen.csv')]
df_pivot_all = pd.DataFrame()
print(len(filelist))
i = 0
for filepath_i in filelist:
    i+=1
    print(f'{i} of {len(filelist)}', filepath_i)
    df = pd.read_csv(filepath_i)
    df['LENGTH'] = df['END'] - df['START']
    df['REVERSED_MAP'] = ['rev' if x else 'for' for x in df['REVERSED_MAPPING']]
    #     df['CLASS'] = df['ALLELE'].str[:7]

    #     df = df[df['CLASS'].isin(['Mamu-B_','Mamu-A1'])]
    df = df[['NAME', 'START', 'END', 'REVERSED_MAP', 'REF_LEN']].drop_duplicates(subset=['NAME', 'REVERSED_MAP'])
    df_m = df.melt(id_vars=['NAME', 'REVERSED_MAP', 'REF_LEN'],
                   value_vars=['START', 'END'],
                   value_name='data', var_name='label',
                   ignore_index=False)

    df_pivot = df_m.pivot_table(values=['data'],
                                index=['NAME', 'REF_LEN'],
                                columns=['REVERSED_MAP', 'label'],
                                aggfunc=np.max,
                                fill_value=np.nan).reset_index()
    df_pivot.columns = ['|'.join(x).strip('|') for x in df_pivot.columns]


    def frag_len(w, x, y, z):
        if math.isnan(x):
            return y - z
        if math.isnan(y):
            return w - x
        return y - x


    df_pivot['MAPPED_LENGTH'] = [frag_len(w, x, y, z) for w, x, y, z in zip(df_pivot['data|for|END'],
                                                                            df_pivot['data|for|START'],
                                                                            df_pivot['data|rev|END'],
                                                                            df_pivot['data|rev|START'])]
    df_pivot_i = df_pivot.dropna()
    df_pivot_i = df_pivot_i[(df_pivot_i['data|for|START'] != 0) & (df_pivot_i['data|rev|START'] != 0)]
    df_pivot_i = df_pivot_i[
        (df_pivot_i['data|for|END'] != df_pivot_i['REF_LEN']) & (df_pivot_i['data|rev|END'] != df_pivot_i['REF_LEN'])]
    df_pivot_i['file'] = os.path.basename(filepath_i)[:6]
    df_pivot_all = pd.concat([df_pivot_all, df_pivot_i], ignore_index=True)
out_path = os.path.join(out_dir, 'aligned_read_lengths.csv')
df_pivot_all.to_csv(out_path, index=False)
