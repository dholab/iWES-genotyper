#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentTypeError
import os
import pandas as pd
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt
import json
from scipy.stats.stats import pearsonr
from scipy import stats

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
    parser.add_argument('-i',
                        '--in_dirs',
                        nargs='+',
                        default=[],
                        help='space delimited list of directories of depth of complete coverage files',
                        required=False)
    parser.add_argument('-s',
                        '--summary_paths',
                        nargs='+',
                        default=[],
                        help='summary path files full paths or relative paths',
                        required=False)
    parser.add_argument('-o', '--out_dir',
                        type=str,
                        help='out put for summary and correlation coeff.json',
                        required=True)
    parser.add_argument('-p', '--species',
                        type=str,
                        default='MCM',
                        help='Species (currently only MCM, Mafa and Mamu allowed)',
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
in_dirs = args.in_dirs
out_dir = args.out_dir
summary_paths = args.summary_paths
species = args.species
in_dirs = [x.strip() for x in in_dirs if len(x.strip()) > 0]
summary_paths = [x.strip() for x in summary_paths if len(x.strip()) > 0]
print('out_dir', out_dir)
print('in_dirs', in_dirs)
print('summary_paths', summary_paths)
print(len(in_dirs))
file_num = 0
df_summary = pd.DataFrame()
max_median = 125
unique_maps_per_allele = 1.5
if species in ['MCM','Mafa']:
    haplotype_group_list = ['Mafa-A1', 'Mafa-B', 'NO_FILTER']
    max_median = 250
    unique_maps_per_allele = 2.5
elif species == 'Mamu':
    haplotype_group_list = ['Mamu-A1', 'Mamu-B', 'NO_FILTER']
    max_median = 125
    unique_maps_per_allele = 1.5
else:
    haplotype_group_list = ['Mamu-A1', 'Mamu-B', 'NO_FILTER']


if (len(in_dirs) + len(summary_paths)) < 1:
    print('No files given')
    exit(1)
if len(in_dirs) > 0:
    for depth_dir in in_dirs:
        file_num += 1
        file_list = os.listdir(depth_dir)
        file_list = [os.path.join(depth_dir, x) for x in file_list if not x.startswith('._') and x.endswith('.csv')]
        df = pd.DataFrame()
        print(len(file_list))
        i = 0
        for file_i in file_list:
            i += 1
            print('file: {0}, file_num: {1} of {2}'.format(file_i, i, len(file_list)))
            df_i = pd.read_csv(file_i)
            df_i = df_i[df_i['HAPLOTYPE_GROUP'].isin(haplotype_group_list)]
            # df_i = df_i[df_i['unique_maps_per_allele'] <1.5]
            df = pd.concat([df, df_i], ignore_index=True)

        df['MAX_POSITION'] = df.groupby(['SAMPLE_NUM', 'ALLELE'])['POSITION'].transform('max')
        df['POSITION_FROM_MAX'] = df['MAX_POSITION'] - df['POSITION']
        df = df[(df['POSITION_FROM_MAX'] >50) & (df['POSITION'] >50)]
        df_summary_i = df.groupby(['SAMPLE_NUM',
                                   'ALLELE',
                                   'HAPLOTYPE_GROUP',
                                   'unique_maps_per_allele']).agg({'DEPTH': ['min', 'median']}).reset_index()

        df_summary_i.columns = ['_'.join(col) for col in df_summary_i.columns.values]
        df_summary_i.columns = [col[:-1] if col.endswith('_') else col for col in df_summary_i.columns]
        # df_summary.columns = [col[:-7] if col.endswith('_median') else col for col in df_summary.columns]

        df_summary_i['SAMPLE_MEDIAN_DEPTH'] = df_summary_i.groupby(['SAMPLE_NUM'])['DEPTH_median'].transform('median')
        df_summary_i['FILE_NUM'] = file_num
        df_summary_i.to_csv(os.path.join(out_dir, f'all_depth_summary_{file_num}.csv'), index=False)
        print('Summary Complete  file_num: {0} of {1}'.format(i, len(file_list)))
        df_summary = pd.concat([df_summary, df_summary_i], ignore_index=True)
for summary_path_i in summary_paths:
    file_num += 1
    df_summary_i = pd.read_csv(summary_path_i)
    df_summary_i['FILE_NUM'] = file_num

    df_summary = pd.concat([df_summary, df_summary_i], ignore_index=True)

df_summary = df_summary[df_summary['SAMPLE_MEDIAN_DEPTH'] < max_median]
df_summary = df_summary[df_summary['unique_maps_per_allele'] < unique_maps_per_allele]
df_summary.to_csv(os.path.join(out_dir, 'summary_correlation.csv'), index=False)
print('Summary Complete')
print('summary exported to {0}'.format(os.path.join(out_dir, 'summary_correlation.csv')))
df_summary['MEDIAN_RATIO'] = df_summary['DEPTH_median'] / df_summary['SAMPLE_MEDIAN_DEPTH']


# remove outlier ratios
df_summary_f = df_summary[(df_summary['MEDIAN_RATIO'] < 1.3) & (df_summary['MEDIAN_RATIO']  >0.7)]
predict_dict = {}
# haplotype_group_list = list(df_summary_f['HAPLOTYPE_GROUP'].unique())
for haplotype_group in haplotype_group_list:
    print(haplotype_group)
    if haplotype_group == 'NO_FILTER':
        df_summary_h = df_summary_f.copy()
    else:
        df_summary_h = df_summary_f[df_summary_f['HAPLOTYPE_GROUP'] == haplotype_group]
    if len(df_summary_h) < 1:
        continue
    lin_regression = LinearRegression().fit(df_summary_h['SAMPLE_MEDIAN_DEPTH'].values.reshape(-1, 1),
                                            df_summary_h['DEPTH_min'].values.reshape(-1, 1))
    # lin_regression

    model_line = lin_regression.predict(df_summary_h['SAMPLE_MEDIAN_DEPTH'].values.reshape(-1, 1))
    # model_line
    plt.scatter(df_summary_h['SAMPLE_MEDIAN_DEPTH'].values.reshape(-1, 1), df_summary_h['DEPTH_min'])
    plt.plot(df_summary_h['SAMPLE_MEDIAN_DEPTH'].values.reshape(-1, 1),model_line)
    plt.xlabel("Median Depth per Sample")
    plt.ylabel("Minimum Depth of Coverage per Allele")
    plt.savefig(os.path.join(out_dir, '{0}_summary_correlation.png'.format(haplotype_group)))
    plt.clf()
    sum_errs = np.sum((df_summary_h['DEPTH_min'].values.reshape(-1, 1)- model_line)**2)
    stdev = np.sqrt(1 / (len(df_summary_h['DEPTH_min'].values.reshape(-1, 1)) - 2) * sum_errs)
    predict_dict[haplotype_group] = [lin_regression.coef_[0][0].round(3),
                                     lin_regression.intercept_[0].round(3),
                                     stdev.round(3)]

    print('CorrelationPlot exported to {0}'.format(os.path.join(out_dir,
                                                                '{0}_summary_correlation.png'.format(haplotype_group))))
# Data to be written
with open(os.path.join(out_dir, "confidence_coeff.json"), "w") as outfile:
    json.dump(predict_dict, outfile)

print('Linear correlation coeff\'s (slope, int, stdev) exported to {0}'.format(os.path.join(out_dir,
                                                                                            'confidence_coeff.json')))
