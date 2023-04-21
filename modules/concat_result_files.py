#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentTypeError
import os
import pandas as pd

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
                        help='directory to output the results',
                        required=True)
    parser.add_argument('--submit_name',
                        type=str,
                        help='Project name for file prefix',
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
in_dir = args.in_dir
submit_name = args.submit_name

df_norm_median = pd.DataFrame()
df_read_ct = pd.DataFrame()
df_gaps = pd.DataFrame()
file_list = os.listdir(in_dir)
file_list = [x for x in file_list if not x.startswith('._') and not x.startswith(submit_name)]
median_all_list = [x for x in file_list if x.endswith('_norm_median_all.csv')]
read_ct_list = [x for x in file_list if x.endswith('_read_ct.csv')]
gaps_list = [x for x in file_list if x.endswith('_gaps.csv')]

for median_i in median_all_list:
    try:
        df_norm_median_i = pd.read_csv(os.path.join(in_dir, median_i))
        df_norm_median = pd.concat([df_norm_median, df_norm_median_i], ignore_index=True)
    except pd.errors.EmptyDataError:
        print('Note: {0} was empty. Skipping.'.format(os.path.join(in_dir, median_i)))
        continue


for read_ct_i in read_ct_list:
    try:
        df_read_ct_i = pd.read_csv(os.path.join(in_dir, read_ct_i))
        df_read_ct = pd.concat([df_read_ct, df_read_ct_i], ignore_index=True)
    except pd.errors.EmptyDataError:
        print('Note: {0} was empty. Skipping.'.format(os.path.join(in_dir, read_ct_i)))
        continue
for gaps_i in gaps_list:
    try:
        df_gaps_i = pd.read_csv(os.path.join(in_dir, gaps_i))
        df_gaps = pd.concat([df_gaps, df_gaps_i], ignore_index=True)
    except pd.errors.EmptyDataError:
        print('Note: {0} was empty. Skipping.'.format(os.path.join(in_dir, gaps_i)))
        continue
df_norm_median.to_csv(os.path.join(in_dir, '{0}_norm_median_all.csv'.format(submit_name)), index=False)
df_read_ct.to_csv(os.path.join(in_dir, '{0}_read_ct.csv'.format(submit_name)), index=False)
df_gaps.to_csv(os.path.join(in_dir, '{0}_gaps.csv'.format(submit_name)), index=False)
