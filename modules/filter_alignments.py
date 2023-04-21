#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentTypeError
import os
import pysam
import pandas as pd
import time
import subprocess
import numpy as np
import re
from io import StringIO
import multiprocessing
from functools import partial
import gc


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')


def str2list(v):
    if isinstance(v, str):
        return v.split(',')
    else:
        raise ArgumentTypeError('Comma delimited str value expected.')


def str2ram(v):
    if isinstance(v, int):
        return v
    if v.lower().endswith in ('gb'):
        return int(v[:-2]) * 1000
    if v.lower().endswith in ('mb'):
        return int(v[:-2])
    if v.lower().endswith in ('tb'):
        return int(v[:-2]) * 1000000
    if v.lower().endswith in ('kb'):
        return int(float(v[:-2]) / 1000)
    else:
        raise ArgumentTypeError('Ram value expected (int in mb, or ends with KB, MB, GB, or TB.')


def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument('--out_dir',
                        type=str,
                        help='directory to output the results',
                        required=True)
    parser.add_argument('--bam_dir',
                        type=str,
                        help='directory your input bam files reside',
                        required=True)
    parser.add_argument('--project_name',
                        type=str,
                        help='Project name for file prefix',
                        required=True)
    parser.add_argument('--bait_fasta',
                        type=str,
                        help='Where you genomic-exon2-miseq combined fasta exists',
                        default=None,
                        required=True)
    parser.add_argument('--ipd_ref_hash',
                        type=str,
                        help='filepath where the ipd reference hash (.csv.gz) files exist',
                        default=None,
                        required=True)

    parser.add_argument('--unpaired_edge_threshold',
                        type=int,
                        help='unpaired_edge_threshold how far (bp) from the edge of reference sequence \
                        to allow the paired read to be unmapped (unmated) and still included',
                        default=500,
                        required=False)
    parser.add_argument('--depth_threshold',
                        type=int,
                        help='depth_threshold how minimum depth of coverage for all \
                        (except the near the edge by defied as edge distance threshold)',
                        default=3,
                        required=False)
    parser.add_argument('--edge_distance_threshold',
                        type=int,
                        help='how far from the edge of the ipd reference sequence to allow for zero depth of coverage \
                             depth of coverage does not apply close to the edge, \
                             short sequence miseq db sequences do not apply as it is set to zero.',
                        default=50,
                        required=False)
    parser.add_argument('--hash_chunks_n',
                        type=int,
                        help='The hashtable processes better in chunks to minimize ram consumption',
                        default=20,
                        required=False
                        )
    parser.add_argument('-e',
                        '--db_ext',
                        nargs='+',
                        default=['gen','exon','miseq'],
                        help='extensions to add at end of the headers to descriminate files (i.e. -e gen exon miseq',
                        required=False)
    parser.add_argument('--edge_ignore_length',
                        type=int,
                        help='The maximum ref seq length to ignore edge filter, instead of doing by data types',
                        default=500,
                        required=False
                        )
    parser.add_argument('--expand_with_hash',
                        type=str2bool,
                        help='How to expand a hash',
                        default=True,
                        required=False)
    parser.add_argument('--threads',
                        type=int,
                        help='many threads to use in multiprocessed operations',
                        default=1,
                        required=False)
    parser.add_argument('--output_depth_all',
                        type=str2bool,
                        help='output a depth of coverage table, even if not meeting filter criteria',
                        default=False,
                        required=False)
    parser.add_argument('--mask_path',
                        type=str,
                        help='Path to mask path leave blank if not in use',
                        default=None,
                        required=False)
    parser.add_argument('--use_mask',
                        type=str,
                        help='Path to mask path leave blank if not in use',
                        default=True,
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


# Get arguments and assign any non-required baths
args = get_process_arrays_args()
# same arguments to a local variable by same name as the argument
out_dir = args.out_dir
bam_dir = args.bam_dir
project_name = args.project_name
unpaired_edge_threshold = args.unpaired_edge_threshold
depth_threshold = args.depth_threshold
edge_distance_threshold = args.edge_distance_threshold
hash_chunks_n = args.hash_chunks_n
bait_fasta = args.bait_fasta
ipd_ref_hash = args.ipd_ref_hash
db_ext = args.db_ext
mask_path = args.mask_path
print(db_ext)
expand_with_hash = args.expand_with_hash
edge_ignore_length = args.edge_ignore_length
threads = args.threads
output_depth_all = args.output_depth_all
use_mask = args.use_mask
os.makedirs(out_dir, exist_ok=True)


def filter_count_map_bam(df_bam_m, allele_list):
    df_bam_m = df_bam_m[df_bam_m['ALLELE'].isin(allele_list)]
    df_bam_m['NAME_COUNTS'] = 1 - ((df_bam_m['PAIRED'] - 1) / 2)
    df_bam_m['num_read_maps'] = df_bam_m.groupby(['NAME'])['NAME_COUNTS'].transform('sum')
    df_bam_m['unique_maps_per_allele'] = df_bam_m.groupby(['ALLELE'])['num_read_maps'].transform('min')
    df_bam_m.sort_values(by=['ALLELE', 'START', 'END'], inplace=True)
    return df_bam_m


def process_df_bam(df_bam, ipd_length_dict, sample_i, unpaired_edge_threshold):
    df_bam['SAMPLE_NUM'] = sample_i
    df_bam["DISTANCE_TO_END"] = [ipd_length_dict[y] - x for x, y in zip(df_bam['END'], df_bam['ALLELE'])]
    df_bam["INCLUDE_UNPAIRED"] = [
        1 if ((x < unpaired_edge_threshold) and z) or ((y < unpaired_edge_threshold) and not (z)) else 0 for x, y, z in
        zip(df_bam["START"], df_bam["DISTANCE_TO_END"], df_bam['REVERSED_MAPPING'])]
    df_bam = df_bam[(df_bam['PAIRED'] + df_bam["INCLUDE_UNPAIRED"]) > 1]
    df_bam['NAME_COUNTS'] = 1 - ((df_bam['PAIRED'] - 1) / 2)
    # df_bam_m['num_read_maps'] = df_bam_m.groupby(['NAME','ALLELE'])['NAME_COUNTS'].transform('sum')
    df_bam['num_read_maps'] = df_bam.groupby(['NAME'])['NAME_COUNTS'].transform('sum')
    df_bam['unique_maps_per_allele'] = df_bam.groupby(['ALLELE'])['num_read_maps'].transform('min')
    return df_bam


def fasta_length_dict(ipd_fasta_path):
    ipd_length_dict = {}

    genome_fasta_open = pysam.Fastafile(ipd_fasta_path)
    for ref in genome_fasta_open.references:
        ipd_length_dict[ref] = genome_fasta_open.get_reference_length(ref)
    genome_fasta_open.close()

    return ipd_length_dict


def get_gap_allele_list(df_bam_m, max_gap=70):
    max_gap = max_gap + 0.5

    df_bam_m_dedup = df_bam_m[
        ['NAME', 'REVERSED_MAPPING', 'START', 'END', 'PAIRED', 'ALLELE', 'SAMPLE_NUM']].drop_duplicates()
    df_bam_m_dedup['START_ALLELE_MIN'] = df_bam_m_dedup.groupby(['ALLELE', 'SAMPLE_NUM'])['START'].transform('min')
    df_bam_m_dedup['END_ALLELE_MAX'] = df_bam_m_dedup.groupby(['ALLELE', 'SAMPLE_NUM'])['END'].transform('max')

    df_bam_paired = df_bam_m_dedup[df_bam_m_dedup['PAIRED'] == 2]
    # df_bam_unpaired = df_bam_m_dedup[df_bam_m_dedup['PAIRED'] == 1]

    df_bam_paired['START_MIN'] = df_bam_paired.groupby(['ALLELE', 'SAMPLE_NUM', 'NAME'])['START'].transform('min')
    df_bam_paired['END_MAX'] = df_bam_paired.groupby(['ALLELE', 'SAMPLE_NUM', 'NAME'])['END'].transform('max')
    df_bam_paired['INSERT_LEN'] = df_bam_paired['END_MAX'] - df_bam_paired['START_MIN']
    df_bam_overlap = df_bam_paired[df_bam_paired['INSERT_LEN'] < 302]
    df_bam_no_overlap = df_bam_paired[df_bam_paired['INSERT_LEN'] > 301]
    df_bam_overlap['START'] = df_bam_overlap['END_MAX'] - 151
    df_bam_overlap['END'] = df_bam_overlap['START_MIN'] + 151
    df_bam_overlap['START'] = [y if x < y else x for x, y in
                               zip(df_bam_overlap['START_MIN'], df_bam_overlap['START_ALLELE_MIN'])]
    df_bam_overlap['END'] = [y if x > y else x for x, y in
                             zip(df_bam_overlap['END_MAX'], df_bam_overlap['END_ALLELE_MAX'])]

    df_bam_overlap = df_bam_overlap[
        ['NAME', 'REVERSED_MAPPING', 'START', 'END', 'PAIRED', 'ALLELE', 'SAMPLE_NUM']].drop_duplicates()

    df_bam_m_dedup = pd.concat([df_bam_m_dedup, df_bam_overlap], ignore_index=True)

    df_start = df_bam_m_dedup[['ALLELE', 'SAMPLE_NUM', 'START']].drop_duplicates()

    df_start.rename(columns={'START': 'POSITION'}, inplace=True)
    df_start.sort_values(['ALLELE', 'SAMPLE_NUM', 'POSITION'], inplace=True)
    df_start['GAP'] = df_start.groupby(['ALLELE', 'SAMPLE_NUM'])['POSITION'].diff()
    df_end = df_bam_m_dedup[['ALLELE', 'SAMPLE_NUM', 'END']].drop_duplicates()
    df_end.rename(columns={'END': 'POSITION'}, inplace=True)
    df_end.sort_values(['ALLELE', 'SAMPLE_NUM', 'POSITION'], inplace=True)
    df_end = df_end[['ALLELE', 'SAMPLE_NUM', 'POSITION']].drop_duplicates()
    df_end['GAP'] = df_end.groupby(['ALLELE', 'SAMPLE_NUM'])['POSITION'].diff()
    df_gap = pd.concat([df_start, df_end], ignore_index=True)
    df_gap.fillna(0, inplace=True)
    df_gap['MAX_GAP'] = df_gap.groupby(['ALLELE', 'SAMPLE_NUM'])['GAP'].transform('max')
    df_gap.sort_values(['ALLELE', 'SAMPLE_NUM', 'POSITION'], inplace=True)
    df_gap_summary = df_gap[df_gap['GAP'] == df_gap['MAX_GAP']]
    df_gap_summary = df_gap_summary.groupby(['ALLELE', 'SAMPLE_NUM', 'MAX_GAP'])['POSITION'].min().reset_index()
    df_gap_summary.rename(columns={'POSITION': 'GAP_POSITION'}, inplace=True)
    df_gap_filtered = df_gap[df_gap['MAX_GAP'] < max_gap + 0.5]

    allele_list = list(df_gap_filtered['ALLELE'].unique())
    #     df_bam_m = df_bam_m[df_bam_m['ALLELE'].isin(allele_list) & df_bam_m['ALLELE_2'].isin(allele_list)]
    return allele_list, df_gap, df_gap_summary


def make_range_from_fasta_lengths_dict(ipd_length_dict, allele_list):
    df_depth_range = pd.DataFrame()
    for allele_i in allele_list:
        df_depth_range_i = pd.DataFrame({'ALLELE': allele_i, 'POSITION': list(range(0, ipd_length_dict[allele_i]))})
        df_depth_range = pd.concat([df_depth_range, df_depth_range_i], ignore_index=True)
    return df_depth_range


def get_depth_and_filter(bam_filepath,
                         edge_distance_threshold=50,
                         depth_threshold=3,
                         sample_i=None,
                         db_ext_i=None,
                         out_dir='./',
                         output_depth_all=False,
                         df_range_mask=pd.DataFrame()):
    # Just use samtools to run the depth. We have a valid Bam file to work off, it is compiled in C
    cmd = ['samtools', 'depth', '-a',
           '-o', '-',
           bam_filepath]
    a = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    b = StringIO(a.communicate()[0].decode('utf-8'))

    df_depth = pd.read_csv(b,
                           sep="\t",
                           header=None,
                           names=['ALLELE', 'POSITION', 'DEPTH'])
    if output_depth_all:
        depth_dir = os.path.join(out_dir, 'depth_all', '{0}.depth_all{1}.csv'.format(sample_i, db_ext_i))
        df_depth.to_csv(depth_dir, sep=',', index=False)
    df_depth['MAX_POSITION'] = df_depth.groupby(['ALLELE'])['POSITION'].transform('max')
    # if the df_range_mask exists then use it to filter
    if len(df_range_mask) > 0:
        allele_list = list(df_range_mask['ALLELE'].unique())
        # figur out the ranges for each allele
        lows_dict = {}
        ups_dict = {}
        for allele_i in allele_list:
            lows_dict[allele_i] = np.array(list(df_range_mask[allele_i == df_range_mask['ALLELE']]['START']))
            ups_dict[allele_i] = np.array(list(df_range_mask[allele_i == df_range_mask['ALLELE']]['END']))

        def in_range(x, lows, ups):
            return np.any((lows <= x) & (x <= ups))

        # calculate if it should be excluded or not from the depth of coverage
        # Exclude for positions means a zero depth of coverage is passes the filter for those positions:
        df_depth['EXCLUDE'] = [in_range(x, lows_dict[y], ups_dict[y]) for x, y in
                               zip(df_depth['POSITION'], df_depth['ALLELE'])]
    else:
        # if there is no range filter set everything to False, meaning, consider every position to be unmasked.
        df_depth['EXCLUDE'] = False
    df_depth = df_depth[(df_depth['EXCLUDE']) |
                        (edge_distance_threshold >= df_depth['POSITION']) |
                        (edge_distance_threshold >= (df_depth['MAX_POSITION'] - df_depth['POSITION'])) |
                        (depth_threshold <= df_depth['DEPTH'])
                        ]
    df_depth['ALLELE_COUNT'] = df_depth.groupby(['ALLELE'])['DEPTH'].transform('count')
    df_depth = df_depth[(df_depth['ALLELE_COUNT'] >= df_depth['MAX_POSITION'])]
    allele_list = list(df_depth['ALLELE'].unique())
    del df_depth
    return allele_list


def depth_worker(allele_chunk):
    df_depth_range_i = allele_chunk['df_depth_range']
    df_bam_m_dedup_i = allele_chunk['df_bam_dedup']
    df_bam_m_dedup_i['POSITION'] = [list(range(x, y)) for x, y in
                                    zip(df_bam_m_dedup_i['START'], df_bam_m_dedup_i['END'])]
    df_depth_i = df_bam_m_dedup_i.explode('POSITION')
    # free up a little bit of ram from variable no longer needed.
    del df_bam_m_dedup_i
    gc.collect()
    df_depth_i = df_depth_i.groupby(['ALLELE', 'POSITION', 'unique_maps_per_allele']).agg(
        {'START': 'count', 'num_read_maps': 'sum'}).rename(columns={'START': 'DEPTH'}).reset_index()
    df_depth_i['DEPTH_ADJ'] = df_depth_i['DEPTH'] * df_depth_i['DEPTH'] / df_depth_i['num_read_maps']
    df_depth_i['DEPTH_ADJ'] = df_depth_i['DEPTH_ADJ'].round(2)
    df_depth_i['DEPTH_RATIO'] = df_depth_i['num_read_maps'] / df_depth_i['DEPTH']
    df_depth_i['DEPTH_RATIO'] = df_depth_i['DEPTH_RATIO'].round(2)

    df_depth_i = df_depth_i.merge(df_depth_range_i, on=['ALLELE', 'POSITION'], how='outer')
    df_depth_i['HAPLOTYPE_GROUP'] = [x.split('_')[0] for x in df_depth_i['ALLELE']]
    gc.collect()
    return df_depth_i


def get_depth(df_bam_m, ipd_length_dict, threads=1):
    sample_i = list(df_bam_m['SAMPLE_NUM'].unique())[0]
    df_bam_m_dedup = df_bam_m[
        ['ALLELE', 'NAME', 'START', 'END', 'unique_maps_per_allele', 'num_read_maps']].drop_duplicates()
    # print(df_bam_m_dedup.memory_usage(deep=True) / 1e6)
    allele_list = list(df_bam_m_dedup['ALLELE'].unique())
    start = 0
    step = 10
    end = len(allele_list)
    # print(start, step, end)
    # df_depth = pd.DataFrame()
    df_depth_range = make_range_from_fasta_lengths_dict(ipd_length_dict, allele_list)
    allele_chunk_list = []
    for i in range(start, end, step):
        # allele_chunk_list.append(allele_list[i:i + step])
        df_depth_range_i = df_depth_range[df_depth_range['ALLELE'].isin(allele_list[i:i + step])]
        df_bam_m_dedup_i = df_bam_m_dedup[df_bam_m_dedup['ALLELE'].isin(allele_list[i:i + step])]
        allele_chunk_list.append({'df_depth_range': df_depth_range_i, 'df_bam_dedup': df_bam_m_dedup_i})


    # for idx, chunk in enumerate(np.array_split(df_hash, hash_chunks_n)):
    #     chunk_list.append(chunk)
    pool = multiprocessing.Pool(processes=threads)
    result = pool.map(depth_worker, allele_chunk_list)
    pool.close()
    pool.terminate()
    df_depth = pd.concat(result, ignore_index=True)
    df_depth.fillna(value=0, inplace=True)
    df_depth['SAMPLE_NUM'] = sample_i
    df_depth['unique_maps_per_allele'] = df_depth.groupby(['ALLELE'])['unique_maps_per_allele'].transform('max')
    return df_depth


def get_normalized_median_by_allele_sample(df_depth, median_groups, ipd_length_dict=[], depth_threshold=1):
    if len(ipd_length_dict) > 0:
        df_depth_copy = df_depth.copy()
        df_depth_copy = df_depth_copy[df_depth_copy['DEPTH'] >= depth_threshold]
        df_summary = df_depth_copy.groupby(['ALLELE', 'SAMPLE_NUM', 'HAPLOTYPE_GROUP']).agg(
            {'unique_maps_per_allele': 'median',
             'DEPTH_ADJ': 'median',
             'POSITION': ['max', 'min']}).reset_index()
        df_summary.columns = ['_'.join(col) for col in df_summary.columns.values]
        df_summary.columns = [col[:-1] if col.endswith('_') else col for col in df_summary.columns]
        df_summary.columns = [col[:-7] if col.endswith('_median') else col for col in df_summary.columns]
        df_summary['Depth_Threshold_From_End_Position'] = [ipd_length_dict[y] - x - 1 for x, y in
                                                           zip(df_summary['POSITION_max'], df_summary['ALLELE'])]
        df_summary['Depth_Threshold_From_Start_Position'] = df_summary['POSITION_min']
        df_summary.drop(columns=['POSITION_max', 'POSITION_min'], inplace=True)
    else:
        df_summary = df_depth.groupby(['ALLELE', 'SAMPLE_NUM', 'HAPLOTYPE_GROUP']).agg(
            {'unique_maps_per_allele': 'median', 'DEPTH_ADJ': 'median'}).reset_index()
    if len(median_groups) > 0:
        df_summary['HAPLOTYPE_GROUP_temp'] = [x.split('-')[1] if len(x.split('-')) > 1 else x for x in
                                              df_summary['HAPLOTYPE_GROUP']]
        df_summary_median = df_summary[df_summary['HAPLOTYPE_GROUP_temp'].isin(median_groups)]
        df_summary_median.drop(columns=['HAPLOTYPE_GROUP_temp'], inplace=True)
    else:
        df_summary_median = df_summary
    median_depth = df_summary_median['DEPTH_ADJ'].median()
    df_summary['DEPTH_NORM'] = df_summary['DEPTH_ADJ'] / median_depth
    df_summary['DEPTH_NORM'] = df_summary['DEPTH_NORM'].round(3)
    return df_summary


def bam_to_df(bam_path):
    samfile = pysam.AlignmentFile(bam_path, 'rb')
    allele_list = []
    name_list = []
    # paired_list = []
    start_list = []
    end_list = []
    reverse_list = []
    # read_len_list = []
    mapping_quality_list = []
    reference_start_list = []
    query_sequence_list = []
    query_qualities_list = []
    reference_end_list = []
    query_flag_list = []
    query_length_list = []
    for segment in samfile:
        if segment.is_unmapped:
            continue
        allele_list.append(segment.reference_name)
        name_list.append(segment.query_name)
        start_list.append(segment.reference_start)
        end_list.append(segment.reference_end)
        # reverse_list.append(segment.is_reverse)
        # read_len_list.append(segment.query_length)
        if not isinstance(segment.reference_end, int):
            reverse_list.append(None)
        else:
            reverse_list.append(segment.is_reverse)
        mapping_quality_list.append(segment.mapping_quality)
        reference_start_list.append(segment.reference_start)
        reference_end_list.append(segment.reference_end)
        query_sequence_list.append(segment.query_sequence)
        query_qualities_list.append(segment.to_string().split('\t')[10])
        query_flag_list.append(segment.to_string().split('\t')[1])
        query_length_list.append(segment.query_length)
    #         if not segment.is_paired or segment.mate_is_unmapped or segment.is_duplicate:
    #             paired_list.append(0)
    #         else:
    #             paired_list.append(1)
    segment_count = len(name_list)
    # 0 qname
    # 1 flag
    # 2 ref name
    # 3 ref start
    # 4 mapping quality
    # 5 cigar string
    # 6 ref id of mate
    # 7 next ref start
    # 8 inferred insert size overlap (reverse is negative, (END2- START1) non paired is 0)
    # 9 query sequence
    # 10 query sequence score
    df_bam_in = pd.DataFrame({'NAME': name_list,
                              'ALLELE_2': allele_list,
                              'START_2': start_list,
                              'END_2': end_list,
                              'REVERSED_MAPPING': reverse_list,
                              'READ_LEN': query_length_list,
                              'MAPPING_QUALITY': mapping_quality_list,
                              'REFERENCE_START': reference_start_list,
                              'REFERENCE_END': reference_end_list,
                              'QUERY_SEQUENCE': query_sequence_list,
                              'QUERY_QUALITIES': query_qualities_list,
                              'QUERY_FLAG': query_flag_list
                              })
    df_bam_in.sort_values(by=['ALLELE_2', 'START_2', 'END_2'], inplace=True)
    return df_bam_in, segment_count


def calc_flag(read_mapped, mate_unmapped=True, reverse=True, first_in_pair=1, read_paired=1):
    # assumes all paired reads
    # read map is 4, mate unmapped is 8 (4+4)
    # reverse is 16/32 (4+16 + 64)
    # 2 is read mapped in proper pair
    # 1	1	Read paired
    # 2	2	Read mapped in proper pair
    # 3	4	Read unmapped
    # 4	8	Mate unmapped
    # 5	16	Read reverse strand
    # 6	32	Mate reverse strand
    # 7	64	First in pair
    # 8	128	Second in pair
    # 9	256	Not primary alignment
    # 10	512	Read fails platform/vendor quality checks
    # 11	1024	Read is PCR or optical duplicate
    # 12	2048	Supplementary alignment
    flag = read_paired
    if read_mapped and not mate_unmapped:
        flag += 2
    if not read_mapped:
        flag += 4
    if mate_unmapped:
        flag += 8
    if reverse:
        flag += 16
    if not mate_unmapped and not reverse:
        flag += 32
    if first_in_pair == 1:
        flag += 64
    if first_in_pair > 1:
        flag += 128
    return flag


def cigar_string(start, end, read_len):
    seg_len = end - start
    if seg_len < read_len:
        split_len = read_len - seg_len
        if start > 0:
            return '{0}={1}S'.format(int(seg_len), int(split_len))
        return '{0}S{1}='.format(int(split_len), int(seg_len))
    return '{0}='.format(read_len)


def makebam(df_bam_in, df_mapped, ipd_length_dict, out_sam):
    df_read_names = df_mapped[['NAME']].drop_duplicates()
    # df_bam_in = df_bam_in[['NAME', 'REVERSED_MAPPING','MAPPING_QUALITY',
    #                        'QUERY_SEQUENCE',
    #                        'QUERY_QUALITIES',
    #                        'READ_LEN']].drop_duplicates()
    df_mapped_dedup = df_mapped[['NAME', 'REVERSED_MAPPING', 'START', 'END', 'PAIRED', 'ALLELE']].drop_duplicates()
    df_mapped_bam = df_bam_in.merge(df_read_names, on=['NAME'], how='inner')
    df_mapped_bam = df_mapped_bam.merge(
        df_mapped_dedup[['NAME', 'REVERSED_MAPPING', 'START', 'END', 'PAIRED', 'ALLELE']],
        on=['NAME', 'REVERSED_MAPPING'],
        how='left')
    df_unmapped = df_mapped_bam[df_mapped_bam['PAIRED'].isna()]
    # Get mapped as they are treated differently
    df_unpaired = df_mapped_bam[df_mapped_bam['PAIRED'] == 1]
    # Get mapped as they are treated differently
    df_paired = df_mapped_bam[df_mapped_bam['PAIRED'] == 2]

    df_unmapped['CIGAR'] = '*'
    df_unmapped['REVERSED_MAPPING'] = False
    df_unmapped['FIRST_IN_PAIR'] = 2
    df_unmapped['READ_MAPPED'] = False
    df_unmapped['MATE_UNMAPPED'] = False
    df_unmapped['MAPPING_QUALITY'] = 0
    df_unmapped['FLAG'] = [calc_flag(read_mapped=w,
                                     mate_unmapped=x,
                                     reverse=y,
                                     first_in_pair=z,
                                     read_paired=1) for w, x, y, z in zip(df_unmapped['READ_MAPPED'],
                                                                          df_unmapped['MATE_UNMAPPED'],
                                                                          df_unmapped['REVERSED_MAPPING'],
                                                                          df_unmapped['FIRST_IN_PAIR'])]

    df_unmapped.drop(columns=['ALLELE', 'START', 'END'], inplace=True)
    df_unpaired_m = df_unpaired[['NAME', 'ALLELE', 'START', 'END']]
    df_unmapped = df_unmapped.merge(df_unpaired_m, on=['NAME'], how='inner')
    df_unmapped['INFERRED_OVERLAP'] = 0

    # df_ipd_unpaired['REVERSED_MAPPING'] = False
    df_unpaired['CIGAR'] = [cigar_string(start, end, read_len) for start, end, read_len in zip(df_unpaired['START'],
                                                                                               df_unpaired['END'],
                                                                                               df_unpaired['READ_LEN'])]
    df_unpaired['FIRST_IN_PAIR'] = 1
    df_unpaired['READ_MAPPED'] = True
    df_unpaired['MATE_UNMAPPED'] = True
    df_unpaired['INFERRED_OVERLAP'] = 0
    df_unpaired['FLAG'] = [calc_flag(read_mapped=w,
                                     mate_unmapped=x,
                                     reverse=y,
                                     first_in_pair=z,
                                     read_paired=1) for w, x, y, z in zip(df_unpaired['READ_MAPPED'],
                                                                          df_unpaired['MATE_UNMAPPED'],
                                                                          df_unpaired['REVERSED_MAPPING'],
                                                                          df_unpaired['FIRST_IN_PAIR'])]

    df_paired['CIGAR'] = [cigar_string(start, end, read_len) for start, end, read_len in zip(df_paired['START'],
                                                                                             df_paired['END'],
                                                                                             df_paired['READ_LEN'])]
    df_paired['FIRST_IN_PAIR'] = df_paired.groupby(['ALLELE', 'NAME'])['START'].rank(method='first')
    df_paired['READ_MAPPED'] = True
    df_paired['MATE_UNMAPPED'] = False
    df_paired['START_1'] = df_paired.groupby(['ALLELE', 'NAME'])['START'].transform('min')
    df_paired['END_1'] = df_paired.groupby(['ALLELE', 'NAME'])['END'].transform('min')
    df_paired['INFERRED_OVERLAP'] = df_paired['END_1'] - df_paired['START_1']
    df_paired['FLAG'] = [calc_flag(read_mapped=w,
                                   mate_unmapped=x,
                                   reverse=y,
                                   first_in_pair=z,
                                   read_paired=1) for w, x, y, z in zip(df_paired['READ_MAPPED'],
                                                                        df_paired['MATE_UNMAPPED'],
                                                                        df_paired['REVERSED_MAPPING'],
                                                                        df_paired['FIRST_IN_PAIR'])]

    df_all_pair = pd.concat([df_paired, df_unpaired, df_unmapped], ignore_index=True)
    allele_list = list(set(df_all_pair['ALLELE']))
    allele_list.sort()

    def concat_headers(allele, length):
        return '@SQ\tSN:{0}\tLN:{1}\n'.format(allele, int(length))

    # zxcv
    # out_bam = '{0}bam'.format(out_sam[:-3])
    print(out_sam)
    with open(out_sam, mode='w') as fout:
        header_string = '@HD\tVN:1.4\tSO:unsorted\n'
        fout.write(header_string)
        # for allele_i in allele_list:
        df_header = pd.DataFrame(ipd_length_dict.items(), columns=['allele', 'length'])
        df_header['length'] = df_header['length'].astype(str)
        df_header = df_header[df_header['allele'].isin(allele_list)]
        # df_header['SAMLINE'] = df_header.apply(lambda x: concat_headers(x['allele'], x['length']), axis=1)
        df_header['SAMLINE'] = '@SQ\tSN:' + df_header['allele'] + '\tLN:' + df_header['length'] + '\n'
        # header_string = header_string + '@SQ\tSN:{0}\tLN:{1}\n'.format(allele_i, int(ipd_length_dict[allele_i]))
        header_string = ''.join(df_header['SAMLINE'])
        fout.write(header_string)
        df_all_pair['START'] = df_all_pair['START'] + 1
        valid_columns = ['NAME',
                         'FLAG',
                         'ALLELE',
                         'START',
                         'CIGAR',
                         'MAPPING_QUALITY',
                         'CIGAR',
                         'INFERRED_OVERLAP',
                         'QUERY_SEQUENCE',
                         'QUERY_QUALITIES']

        dtypes = df_all_pair[valid_columns].dtypes.to_dict()
        for col_name, typ in dtypes.items():
            if (typ == 'int64'):
                df_all_pair[col_name] = df_all_pair[col_name].astype(str)
            if (typ == 'float64'):
                df_all_pair[col_name] = df_all_pair[col_name].astype(int)
                df_all_pair[col_name] = df_all_pair[col_name].astype(str)
        df_all_pair['SAMLINE'] = df_all_pair['NAME'] + '\t' + \
                                 df_all_pair['FLAG'] + '\t' + \
                                 df_all_pair['ALLELE'] + '\t' + \
                                 df_all_pair['START'] + '\t' + \
                                 df_all_pair['MAPPING_QUALITY'] + '\t' + \
                                 df_all_pair['CIGAR'] + '\t' + \
                                 df_all_pair['ALLELE'] + '\t' + \
                                 df_all_pair['START'] + '\t' + \
                                 df_all_pair['INFERRED_OVERLAP'] + '\t' + \
                                 df_all_pair['QUERY_SEQUENCE'] + '\t' + \
                                 df_all_pair['QUERY_QUALITIES'] + '\n'

        sam_lines = ''.join(df_all_pair['SAMLINE'])
        fout.write(sam_lines)
    # 0 qname
    # 1 flag
    # 2 ref name
    # 3 ref start
    # 4 mapping quality
    # 5 cigar string
    # 6 ref id of mate
    # 7 next ref start
    # 8 inferred insert size overlap (reverse is negative, (END2- START1) non paired is 0)
    # 9 query sequence
    # 10 query sequence score

    out_bam = '{0}bam'.format(out_sam[:-3])
    create_bam_cmd = ['samtools', 'view', '-o', '{0}.bam'.format(out_bam), out_sam]
    subprocess.run(create_bam_cmd)

    sort_bam_cmd = ['samtools', 'sort',
                    '-o', out_bam, '{0}.bam'.format(out_bam)]
    subprocess.run(sort_bam_cmd)
    # sort_bam_cmd = 'samtools sort -o {0} {0}.bam'.format(out_bam)
    # os.system(sort_bam_cmd)
    os.remove(out_sam)
    os.remove('{0}.bam'.format(out_bam))
    return out_bam


def annotate_depth(x, y, z):
    y = str(int(y))
    y = y.zfill(2)
    ambig = '00'
    if z < 0.75:
        ambig = '02'
    return '{0}.{1}{2}'.format(int(x), y, ambig)


def expand_bam_hash_worker(chunk, df_bam_in):
    df_exp_i = df_bam_in.merge(chunk, on=['ALLELE_2'], how='inner')
    # save memory by deleting variables not in use
    del chunk
    gc.collect()
    # and replace negative with 0
    # 2_x is bam in
    # 2-Y is look up on allele2
    df_exp_i['START_NEW'] = (df_exp_i['START_2_x'] - df_exp_i['START_2_y'] + df_exp_i['START']).clip(0)
    df_exp_i['END_NEW'] = df_exp_i['END_2_x'] - df_exp_i['START_2_y'] + df_exp_i['START']
    # Fast method for:    df_exp_i['START_NEW'] = [0 if x < 0 else x for x in  df_exp_i['START_NEW']]
    df_exp_i['END_NEW'] = df_exp_i['END_NEW'] - (df_exp_i['END_NEW'] - df_exp_i['REF_LEN_2']).clip(0)
    df_exp_i['END_NEW'] = df_exp_i['END_NEW'] - (df_exp_i['END_NEW'] - df_exp_i['REF_LEN']).clip(0)
    # The Matching Segment(NEW) must encompass the entire sequence OR
    # The 1st Aligned aligned segment must be within the start and end positions AND
    # ( The matching New segments length must equal the read length.
    # OR both start positions of the New and 1st aligned segment must be zero OR
    # the end positions of the New and 1st aligned segment mus match the  corresponding reference sequence length).

    df_exp_i = df_exp_i[((df_exp_i['START_NEW'] == 0) & (df_exp_i['END_NEW'] == df_exp_i['REF_LEN'])) |
                        ((df_exp_i['END_2_x'] <= df_exp_i['END_2_y']) &
                         (df_exp_i['START_2_x'] >= df_exp_i['START_2_y']) &
                         (df_exp_i['END_NEW'] - df_exp_i['START_NEW'] >= df_exp_i['READ_LEN'])) |
                        (
                         ((df_exp_i['START_NEW'] <= 0) &
                          (df_exp_i['START'] == 0) &
                          (df_exp_i['END_2_x'] <= df_exp_i['END_2_y'])) |
                         ((df_exp_i['END_NEW'] >= df_exp_i['REF_LEN']) &
                          (df_exp_i['END'] == df_exp_i['REF_LEN']) &
                          (df_exp_i['START_2_x'] >= df_exp_i['START_2_y']))
                        )]
    # Also filter out short new match segments less than 1/2 the read segment per semi-perfect matching requirement.
    df_exp_i = df_exp_i[(df_exp_i['READ_LEN'] / 2) < (df_exp_i['END_NEW'] - df_exp_i['START_NEW'])]
    df_exp_i.drop(columns=['START_2_y', 'END_2_y', 'REF_LEN_2', 'START', 'END', 'READ_LEN'], inplace=True)
    df_exp_i.rename(columns={'START_2_x': 'START_2', 'END_2_x': 'END_2', 'START_NEW': 'START', 'END_NEW': 'END'},
                    inplace=True)
    df_exp_i['START_2'] = df_exp_i['START_2'].astype(int)
    df_exp_i['END_2'] = df_exp_i['END_2'].astype(int)
    # The pool probably locks some automatic garbage collection to free memory
    # Also setting a filter onto itself probably needs to create a new object in the background
    gc.collect()

    return df_exp_i


def expand_bam_hash(df_hash, hash_chunks_n, df_bam_in, threads=1):
    chunk_list = []
    pool = multiprocessing.Pool(processes=threads)
    for idx, chunk in enumerate(np.array_split(df_hash, hash_chunks_n)):
        chunk_list.append(chunk)
    result = pool.map(partial(expand_bam_hash_worker, df_bam_in=df_bam_in), chunk_list)
    pool.close()
    pool.terminate()
    df_all = pd.concat(result, ignore_index=True)
    df_all['SEGMENT_LEN'] = df_all['END'] - df_all['START']
    df_all.sort_values(by=['SEGMENT_LEN'], ascending=False, inplace=True)
    df_all.drop_duplicates(subset=['NAME', 'ALLELE_2', 'REVERSED_MAPPING', 'ALLELE'], keep='first', inplace=True)
    df_all.drop(columns=['SEGMENT_LEN'], inplace=True)
    return df_all


# Removes allows us to combine classes in a single command
pd.options.mode.chained_assignment = None
# Get timers to track how long things take
t1 = time.time()
t2 = time.time()

ipd_length_dict = fasta_length_dict(bait_fasta)
# use a hash instead
i = 1
bam_filelist = os.listdir(bam_dir)
bam_filelist = [os.path.join(bam_dir, x) for x in bam_filelist if x.endswith('.bam') and not x.startswith('._')]
df_median_dict = {}
df_gap_dict = {}
# convert to - notation.
db_ext = ['-{0}'.format(x) for x in db_ext]
for db_ext_i in db_ext:
    df_gap_dict[db_ext_i] = pd.DataFrame()
    df_median_dict[db_ext_i] = pd.DataFrame()
df_gap_summary = pd.DataFrame()
# df_depth_flagged = pd.DataFrame()  # Not sure
df_read_ct_all = pd.DataFrame()
bam_count = len(bam_filelist)
bam_iter = 1

# how many chunks for the hash to save memory
# open the files as needed and piece together into nice files
df_hash = pd.DataFrame()

if os.path.isfile(ipd_ref_hash):
    df_hash = pd.read_csv(ipd_ref_hash)
    df_length = pd.DataFrame(list(ipd_length_dict.items()), columns=['ALLELE', 'REF_LEN'])
    df_hash = df_hash.merge(df_length, on=['ALLELE'])
    # print(df_hash)
    df_hash['REF_LEN_2'] = [ipd_length_dict[x] for x in df_hash['ALLELE_2']]
    # print('df_hash', df_hash)

for bam_file_i in bam_filelist:
    print('--- {0} of {1} ---'.format(bam_iter, bam_count))
    bam_iter += 1
    sample_i = os.path.basename(bam_file_i).split('.')[0]
    df_bam_in, segment_count = bam_to_df(bam_file_i)
    df_read_ct_i = pd.DataFrame({'sample_read_ct': [segment_count], 'gs_id': [sample_i]})
    df_read_ct_all = pd.concat([df_read_ct_all, df_read_ct_i], ignore_index=True)
    if len(df_bam_in) < 1:
        print('Bam file could not be converted to a dataframe (could be empty): {0}'.format(sample_i))
        continue
    i = 0
    print(time.time() - t2)
    # hash_chunks_n = 4
    df_bam_in_t = df_bam_in[['NAME', 'ALLELE_2', 'START_2', 'END_2', 'REVERSED_MAPPING', 'READ_LEN']]
    df_all = expand_bam_hash(df_hash, hash_chunks_n, df_bam_in_t, threads)

    del df_bam_in_t
    df_all.drop_duplicates(inplace=True)
    print(time.time() - t2)
    print('length of expanded dt: {0}'.format(len(df_all)))
    print(df_all.dtypes.to_dict())
    df_all['PAIRED'] = df_all.groupby(['NAME', 'ALLELE'])['REVERSED_MAPPING'].transform('count')
    # df_all['START'] = df_all['START'].astype(int)
    # df_all['END'] = df_all['END'].astype(int)
    if use_mask and mask_path is not None and os.path.isfile(mask_path):
        try:
            df_range_mask = pd.read_csv(mask_path)
        except:
            df_range_mask = pd.DataFrame()
    else:
        df_range_mask = pd.DataFrame()
    print('Ensure proper types', time.time() - t2)
    for db_ext_i in db_ext:
        # df_dict[db_ext_i] = df_all[df_all['ALLELE'].str.endswith(db_ext_i)]
        df_ext = df_all[df_all['ALLELE'].str.endswith(db_ext_i)]
        print('Subset grouping', time.time() - t2)
        df_bam = process_df_bam(df_bam=df_ext,
                                ipd_length_dict=ipd_length_dict,
                                sample_i=sample_i,
                                unpaired_edge_threshold=unpaired_edge_threshold)
        del df_ext
        print('Process BAMS', time.time() - t2)
        # Output Raw Bam file that is expanded as needed (Have a wrapper for exceptions i.e. blank file)
        # Remove Unpaired reads that are paired with another viable allele
        df_bam['MAX_PAIRS'] = df_bam.groupby(['NAME'])['PAIRED'].transform('max')
        df_bam = df_bam[df_bam['MAX_PAIRS'] == df_bam['PAIRED']]
        print('Filter for best pairings', time.time() - t2)
        if not ((len(df_bam_in) > 0) and (len(df_bam) > 0)):
            continue
        print('Export bam file including alleles without full coverage: {0}'.format(sample_i))
        os.makedirs(os.path.join(out_dir, 'bam_all{0}'.format(db_ext_i)), exist_ok=True)
        out_bam = makebam(df_bam_in=df_bam_in,
                          df_mapped=df_bam,
                          ipd_length_dict=ipd_length_dict,
                          out_sam=os.path.join(out_dir,
                                               'bam_all{0}'.format(db_ext_i),
                                               '{0}_all{1}.sam'.format(sample_i, db_ext_i)))
        print('Make full covered bam', time.time() - t2)
        print('Finished IPD pair_reads_from_bam: {0}'.format(sample_i))

        # aggregating for depth of coverage is time consuming. We can eliminate ~95% looking for <1 or more. Then we can
        # calculate the depth of coverage for a 20x speed boost.
        allele_list = get_depth_and_filter(bam_filepath=out_bam,
                                           edge_distance_threshold=edge_distance_threshold,
                                           depth_threshold=depth_threshold,
                                           sample_i=sample_i,
                                           db_ext_i=db_ext_i,
                                           out_dir=out_dir,
                                           output_depth_all=output_depth_all,
                                           df_range_mask=df_range_mask)
        print('Filter for depth and make depth coverage using samtools', time.time() - t2)
        df_bam = filter_count_map_bam(df_bam, allele_list)
        print('filter_count_map_bam', time.time() - t2)
        if len(df_bam) < 1:
            print('1st  Depth of coverage (DOC < 0) filter filtered out all alleles for sample: {0}'.format(sample_i))
            continue
        print('Finished quick screening depth 1 of 4: {0}'.format(sample_i))
        df_depth = get_depth(df_bam, ipd_length_dict, threads)
        print('Get depth and normalized depth', time.time() - t2)
        allele_list_, df_gap, df_gap_summary_i = get_gap_allele_list(df_bam)
        print('Calculate gaps', time.time() - t2)
        # df_gap_summary = pd.concat([df_gap_summary, df_gap], ignore_index=True)
        df_gap_dict[db_ext_i] = pd.concat([df_gap_dict[db_ext_i], df_gap], ignore_index=True)
        del df_gap
        os.makedirs(os.path.join(out_dir, 'expanded_maps{0}'.format(db_ext_i)), exist_ok=True)
        df_bam.to_csv(os.path.join(out_dir,
                                   'expanded_maps{0}'.format(db_ext_i),
                                   '{0}_expanded_maps{1}2.csv'.format(sample_i, db_ext_i)),
                      index=False)
        df_bam_2 = df_bam.copy()
        df_bam_2.rename(columns={'ALLELE': 'ALLELE_2'}, inplace=True)
        df_bam_2 = df_bam_2[['NAME', 'REVERSED_MAPPING', 'ALLELE_2']]
        df_bam_m = df_bam.merge(df_bam_2, on=['NAME', 'REVERSED_MAPPING'])
        del df_bam_2
        os.makedirs(os.path.join(out_dir, 'expanded_maps{0}'.format(db_ext_i)), exist_ok=True)
        df_bam_m.to_csv(os.path.join(out_dir,
                                     'expanded_maps{0}'.format(db_ext_i),
                                     '{0}_expanded_maps{1}.csv'.format(sample_i, db_ext_i)),
                        index=False)
        del df_bam_m
        os.makedirs(os.path.join(out_dir, 'depth{0}'.format(db_ext_i)), exist_ok=True)
        df_depth.to_csv(os.path.join(out_dir,
                                     'depth{0}'.format(db_ext_i),
                                     '{0}_depth{1}.csv'.format(sample_i, db_ext_i)), index=False)
        print('Make intermediate files', time.time() - t2)
        df_norm_median_i = get_normalized_median_by_allele_sample(df_depth,
                                                                  median_groups=['A1', 'B'],
                                                                  ipd_length_dict=ipd_length_dict,
                                                                  depth_threshold=depth_threshold)
        del df_depth
        df_norm_median_i = df_norm_median_i.merge(df_gap_summary_i, on=['ALLELE',
                                                                        'SAMPLE_NUM'], how='inner')
        del df_gap_summary_i
        os.makedirs(os.path.join(out_dir, 'normalized_median{0}'.format(db_ext_i)), exist_ok=True)
        df_norm_median_i.to_csv(os.path.join(out_dir,
                                             'normalized_median{0}'.format(db_ext_i),
                                             '{0}_norm_median{1}.csv'.format(sample_i, db_ext_i)),
                                index=False)
        print('Make normalized median files', time.time() - t2)
        print('--- start make filtered and expanded bam files ---')
        if (len(df_bam_in) > 0) and (len(df_bam) > 0):
            os.makedirs(os.path.join(out_dir, 'bam{0}'.format(db_ext_i)), exist_ok=True)

            makebam(df_bam_in=df_bam_in,
                    df_mapped=df_bam,
                    ipd_length_dict=ipd_length_dict,
                    out_sam=os.path.join(out_dir, 'bam{0}'.format(db_ext_i), '{0}{1}.sam'.format(sample_i, db_ext_i)))
        print('--- {0} complete for {1} ---'.format(db_ext_i, sample_i))
        del df_bam

        df_median_dict[db_ext_i] = pd.concat([df_median_dict[db_ext_i], df_norm_median_i], ignore_index=True)
        print('Concat files as needed', time.time() - t2)
        del df_norm_median_i
    del df_bam_in
    print('--- SAMPLE complete {0} ---'.format(sample_i))
    print(time.time() - t2)
    t2 = time.time()
df_norm_median_all = pd.DataFrame()
df_gap_all = pd.DataFrame()
for db_ext_i in db_ext:
    if len(df_median_dict[db_ext_i]) > 0:
        df_median_dict[db_ext_i]['DEPTH_ADJ'] = df_median_dict[db_ext_i]['DEPTH_ADJ'].astype(int)
        df_median_dict[db_ext_i]['DB'] = re.sub(r'\W+', '', db_ext_i)
        df_norm_median_all = pd.concat([df_norm_median_all, df_median_dict[db_ext_i]], ignore_index=True)
        df_gap_dict[db_ext_i]['DB'] = re.sub(r'\W+', '', db_ext_i)
        df_gap_all = pd.concat([df_gap_all, df_gap_dict[db_ext_i]], ignore_index=True)

df_norm_median_all.to_csv(os.path.join(out_dir, '{0}_norm_median_all.csv'.format(project_name)), index=False)
df_read_ct_all.to_csv(os.path.join(out_dir, '{0}_read_ct.csv'.format(project_name)), index=False)
df_gap_all.to_csv(os.path.join(out_dir, '{0}_gap_summary.csv'.format(project_name)), index=False)
# df_depth_flagged.to_csv(os.path.join(out_dir, '{0}_depth_flagged.csv'.format(project_name)), index=False)
print(out_dir)
print(time.time() - t1)
print('---Pipeline Complete--')
