from matplotlib.patches import Patch
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from os.path import join as join_path
from argparse import ArgumentParser, ArgumentTypeError
import pandas as pd
import numpy as np
import os
import math
from sklearn.preprocessing import LabelEncoder
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
                        help='dir path to the expanded_maps produced by this pipeline',
                        required=True)
    parser.add_argument('--out_dir',
                        type=str,
                        help='path to your output dir.',
                        required=True)
    parser.add_argument('--prefix',
                        type=str,
                        default='104294',
                        help='sample name that is prefix for file name of expanded maps such as: 104294',
                        required=False)
    parser.add_argument('-a',
                        '--allele_list',
                        nargs=2,
                        default=['Mamu-A1_028_03_01_01-gen', 'Mamu-A1_028_01_01_01-gen'],
                        help='extensions to add at end of the headers to descriminate files (i.e. -e gen exon miseq',
                        required=False)
    parser.add_argument('-s',
                        '--start',
                        help='extensions to add at end of the headers to descriminate files (i.e. -e gen exon miseq',
                        type=int,
                        required=True)
    parser.add_argument('-e',
                        '--end',
                        type=int,
                        help='extensions to add at end of the headers to descriminate files (i.e. -e gen exon miseq',
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

prefix = args.prefix
out_dir = args.out_dir
in_dir = args.in_dir
start = args.start
end = args.end
allele_list = args.allele_list
le = LabelEncoder()
print(allele_list)
filelist = os.listdir(in_dir)
filelist = [os.path.join(in_dir, x) for x in filelist if
            not x.startswith('._') and x.endswith('_expanded_maps-gen.csv')]
filelist = [x for x in filelist if prefix in x]
df_pivot_all = pd.DataFrame()
print(len(filelist))


def color(row, c_dict):
    return c_dict[row]


for file_i in filelist:
    print(file_i)
    df = pd.read_csv(file_i)
    print(df.columns)
    df['LENGTH'] = df['END'] - df['START']
    df['REVERSED_MAP'] = ['rev' if x else 'for' for x in df['REVERSED_MAPPING']]
    df['TOKEN_NAME'] = le.fit_transform(df['NAME'])
    df['TOKEN_NAME'] = df['TOKEN_NAME'].astype(str)
    print(df.columns)
    df_i = df[['TOKEN_NAME', 'REVERSED_MAP', 'ALLELE', 'REF_LEN', 'START', 'END', 'LENGTH', 'ALLELE_2_y.1']]
    df_i['MAX_END'] = df_i.groupby(['TOKEN_NAME', 'ALLELE', 'ALLELE_2_y.1'])['END'].transform(max)
    df_i['MIN_START'] = df_i.groupby(['TOKEN_NAME', 'ALLELE', 'ALLELE_2_y.1'])['START'].transform(min)
    df_i['MAX_LEN'] = df_i['MAX_END'] - df_i['MIN_START']
    df_i['SEGMENT_LEN'] = [y if y < 286 else x for x, y in zip(df_i['LENGTH'], df_i['MAX_LEN'])]

    # df_i = df_i[~((df_i['REVERSED_MAP']=='rev')& (df_i['SEGMENT_LEN']>=df_i['LENGTH']))]
    df_i['TOKEN_NAME_DIR'] = [f'{x}_{y}' for x, y in zip(df_i['TOKEN_NAME'], df_i['REVERSED_MAP'])]

    df_i = df_i[df_i['ALLELE'].isin(allele_list)]
    df_i = df_i[(df_i['START'] > start) & (df_i['START'] < end)]
    df_i = df_i.drop_duplicates(subset=['TOKEN_NAME_DIR', 'ALLELE', 'REVERSED_MAP', 'ALLELE_2_y.1'])

    df_i['ALLELE_MAPPING_LIST'] = df_i.groupby(['TOKEN_NAME', 'ALLELE', 'REVERSED_MAP'])['ALLELE_2_y.1'].transform(
        lambda x: '\n'.join(x))

    df_i = df_i.drop_duplicates(subset=['TOKEN_NAME', 'ALLELE', 'REVERSED_MAP'])

    allele_mapping_list = list(df_i.ALLELE_MAPPING_LIST.unique())
    print(allele_mapping_list)
    color_j = 0
    color_list = ['#e49e00', '#56b5e9', '#019f73', '#0272b2', '#f0e643']
    c_dict = {}

    for allele_i, color_i in zip(allele_mapping_list, color_list[:len(allele_mapping_list)]):
        c_dict[allele_i] = color_i
    print(c_dict)
    df_i['color'] = df_i.apply(lambda x: color(x['ALLELE_MAPPING_LIST'], c_dict), axis=1)
    for allele_i in allele_list:
        df_j = df_i[(df_i['ALLELE'] == allele_i)]
        df_j.sort_values(by=['START'], ascending=False, inplace=True)
        fig, ax = plt.subplots(1, figsize=(16, 6))
        ax.barh(df_j.TOKEN_NAME_DIR, df_j.LENGTH, left=df_j.START, color=df_j.color)
        legend_elements = [Patch(facecolor=c_dict[i], label=i) for i in c_dict]
        # Hide Y axes tick marks and labels
        ax.yaxis.set_tick_params(labelleft=False)
        ax.set_yticks([])
        plt.xlabel("Position")
        plt.ylabel("Reads")
        plt.legend(handles=legend_elements, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        plt.gcf().set_size_inches(10, 5)
        plt.savefig(os.path.join(out_dir, f'{allele_i}_chimera_example.svg'), dpi=200)
