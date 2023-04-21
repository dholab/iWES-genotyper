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
                        help='where the config files are stored',
                        required=True)
    parser.add_argument('-d',
                        '--db_path',
                        nargs='+',
                        help='Database file paths space separated to add multiple -d /path/one /path/two /path/three',
                        required=True)
    parser.add_argument('-e',
                        '--db_ext',
                        nargs='+',
                        default=None,
                        help='extensions to add at end of the headers to descriminate files (i.e. -e gen exon miseq',
                        required=False)
    parser.add_argument('--haplo_fasta',
                        type=str,
                        help='the database fasta file the haplotype_json is based on',
                        default=None,
                        required=False)
    parser.add_argument('--species',
                        type=str,
                        help='Species prefix',
                        required=True)
    parser.add_argument('--haplotype_json_path',
                        type=str,
                        help='haplotype_json_path',
                        default=None,
                        required=False)
    parser.add_argument('--cp_path',
                        type=str,
                        help='Directory Where your bbmap executables are stored',
                        required=True)
    parser.add_argument('--threads',
                        type=int,
                        help='Number of threads to run bbmap',
                        default=1,
                        required=False)
    parser.add_argument('--ram',
                        type=int,
                        help='Directory ram to dedicate to bbmap java',
                        default=8000,
                        required=False)
    parser.add_argument('--minimap2_path',
                        type=str,
                        help='location of minimap2 executable often it is in ./bin and can be called with minimap2',
                        default='minimap2',
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
ref_dir = args.ref_dir

haplotype_json_path = args.haplotype_json_path
haplo_fasta = args.haplo_fasta
cp_path = args.cp_path
threads = args.threads
ram = args.ram
db_path = args.db_path
db_ext = args.db_ext
print(db_path, db_ext)
species = args.species
minimap2_path = args.minimap2_path
print(minimap2_path)
# change to where your bbmap is to launch java version. (add '/current' to the directory bbmap executables are located)
# what delimiter to split the header names to shorten them: '' is no split
split_char_ipd = ' '
split_char_diag = ''

os.makedirs(ref_dir, exist_ok=True)

# ipd_fasta_dir = os.path.join(ref_dir, 'ipd')
# diag_fasta_dir = os.path.join(ref_dir, 'miseq')
all_filepath = os.path.join(ref_dir, 'bait.fasta')
haplotype_csv_path = os.path.join(ref_dir, 'haplotype_lookup.csv')
pd.options.mode.chained_assignment = None

def rev_comp_st(seq):
    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    seq = seq[::-1]
    return seq


def trim_read_for_diag(seq, cigartupple):
    i = 0
    trim_start = 0
    trim_end = 0
    final_tupple = len(cigartupple) - 1
    for cigar_i in cigartupple:
        if i == 0:
            if cigar_i[0] == 4:
                trim_start = cigar_i[1]
        if i == final_tupple:
            if cigar_i[0] == 4:
                trim_end = cigar_i[1]
        i += 1

    if trim_end == 0:
        return seq[trim_start:]
    return seq[trim_start:-trim_end]


def write_sequence(name, seq, out_file, replace_illegal_char=True, prefix='', suffix=''):
    if replace_illegal_char:
        name = name.replace('*', '_')  # add a double because i want to fine the 02 at the beginning easir
        name = name.replace(':', '_')
    out_file.write(
        '>{0}{1}{2}\n'.format(prefix, name, suffix))  # . We do mostly non human primates I want to differentiate.
    out_file.write('{0}\n'.format(seq))


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


def single_file_per_fasta(fasta_path, single_fasta_dir, name_num_json_path):
    i = 0
    ipd_num = {}
    fasta_sequences = SeqIO.parse(open(fasta_path), 'fasta')
    for fasta in fasta_sequences:
        i += 1
        name, sequence = fasta.id, str(fasta.seq)
        if len(sequence) > 1:
            truncated_ref = os.path.join(single_fasta_dir, '{0}.fasta'.format(i))
            ipd_num[name] = '{0}'.format(i)
            with open(truncated_ref, 'w') as out_file:
                out_file.write('>{0}\n'.format(name))
                out_file.write('{0}\n'.format(sequence))
    with open(name_num_json_path, 'w') as convert_file:
        convert_file.write(json.dumps(ipd_num))
    return


def like_join_haplo(x, df_ref, column_i='SEQUENCE', name_i='SEQUENCE'):
    name_list = []
    for idx, row in df_ref.iterrows():
        #if x in row[column_i]:
        if (x in row[column_i]) or (rev_comp_st(x) in row[column_i]):
            name_list.append(row[name_i])
    return name_list


def like_join(x, df_ref, column_i='SEQUENCE'):
    name_list = []
    for idx, row in df_ref.iterrows():
        if x in row['allele']:
            name_list.append(row[column_i])
    return name_list


def rev_comp_like_join(inner_seq, df_ref):
    name_list = []
    for idx, row in df_ref.iterrows():
        if (inner_seq in row['SEQUENCE']) or (rev_comp_st(inner_seq) in row['SEQUENCE']):
            name_list.append(row['allele'])
    return name_list


#####################################
#   Generate legal name for fasta   #
#####################################
fasta_path_dict = {}
for x, y in zip(db_ext, db_path):
    fasta_path_dict[x] = y
print(fasta_path_dict)
for suffix, path in fasta_path_dict.items():
    fasta_dict = {}
    fasta_sequences = SeqIO.parse(open(path), 'fasta')
    for fasta in fasta_sequences:
        name = fasta.name
        #         if len(name) > name_lengths:
        #             name_lengths = len(name)
        #             name_long = name
        if len(name) > 1:
            fasta_dict[name.split('|')[0]] = str(fasta.seq)
    name_list = list(fasta_dict.keys())
    print(len(name_list), suffix, path)
    name_list.sort()
    with open(os.path.join(ref_dir, '{0}.fasta'.format(suffix)), 'w') as out_file:
        for name in name_list:
            write_sequence(name=name,
                           seq=fasta_dict[name],
                           out_file=out_file,
                           replace_illegal_char=True,
                           prefix='',
                           suffix='-{0}'.format(suffix))

all_paths = [os.path.join(ref_dir, '{0}.fasta'.format(db_ext[0]))]
# concatenate the exon2 and genomic files
for i in range(1, len(db_ext)):
    concat_suffix_list = db_ext[0:i]
    concat_suffix_str = '_'.join(concat_suffix_list)
    # fasta_name = '{0}.fasta'.format(concat_suffix_str)
    concat_fasta_path = os.path.join(ref_dir, '{0}_concat.fasta'.format(concat_suffix_str))
    with open(concat_fasta_path, 'w') as out_file:
        for suffix in concat_suffix_list:
            fasta_sequences = SeqIO.parse(open(os.path.join(ref_dir, '{0}.fasta'.format(suffix))), 'fasta')
            for fasta in fasta_sequences:
                name = fasta.name
                write_sequence(name=name,
                               seq=str(fasta.seq),
                               out_file=out_file,
                               replace_illegal_char=True,
                               prefix='',
                               suffix='')
    compare_fasta_path = fasta_path_dict[db_ext[i]]
    out_bam = os.path.join(ref_dir, f'{db_ext[i]}_to_{concat_suffix_str}.bam')
    # minimap_path = ''
    minimap_cmd = f'{minimap2_path} -N 1 --sam-hit-only --secondary=no -a {compare_fasta_path} {concat_fasta_path} > {out_bam}'
    print(minimap_cmd)
    os.system(minimap_cmd)

    bf_in = pysam.AlignmentFile(out_bam, "r")
    read_list_mapped = []
    # read_list_unmapped = []
    cigar_list_mapped = []
    # sequence_list_map = []
    trimmed_sequence = []
    read_dict = {}
    for r in bf_in.fetch(until_eof=True):
        mapped_read = str(r).split('\t')
        if int(mapped_read[1]) == 0:
            read_list_mapped.append(mapped_read[0])
            cigar_list_mapped.append(r.cigartuples)
            # sequence_list_map.append(r.query_sequence)
            trimmed_sequence.append(r.query_alignment_sequence)
    bf_in.close()
    # sort by name
    df_diag_merge = pd.DataFrame()
    if len(trimmed_sequence) > 0:
        read_list_mapped, trimmed_sequence = zip(*sorted(zip(read_list_mapped, trimmed_sequence)))
        print(trimmed_sequence[:3])
        print(read_list_mapped[:3])

        with open(os.path.join(ref_dir, '{0}_ungrouped.fasta'.format(db_ext[i])), 'w') as out_file:
            for read_name, seq in zip(read_list_mapped, trimmed_sequence):
                write_sequence(name=read_name,
                               seq=seq,
                               out_file=out_file,
                               replace_illegal_char=True,
                               prefix='')

        ###########################################################
        # deduplicate miseq database fasta and name group headers #
        ###########################################################
        print('deduplicate miseq database fasta')
        df_diag = fasta_to_df(fasta_path=os.path.join(ref_dir, '{0}_ungrouped.fasta'.format(db_ext[i])),
                              header_name='allele', sequence_name='SEQUENCE',
                              as_df=True)
        # remove any duplicate sequences as we have them grouped already
        df_diag_dedup = df_diag.drop_duplicates(subset=['SEQUENCE'], keep='first', ignore_index=False)
        df_diag_dedup.rename(columns={'allele': 'allele_first'}, inplace=True)
        # merge to the sequence to find a better grouping name
        # This part renames the grouping but is not very dynamic relies on a lot of convention that may not exist
        df_diag_merge = df_diag.merge(df_diag_dedup, on=['SEQUENCE'], how='inner')
        df_diag_merge['allele'] = ['-'.join(x.split('-')[1:]) if x.count('-') > 1 else x for x in
                                   df_diag_merge['allele']]
        # Turn duplicates to a comma delimited list with aggregation group by.
        # This creates a three column list: "allele_first" = uid, allele = comma list, sequence = sequence
        df_diag_merge = (df_diag_merge.groupby(['SEQUENCE', 'allele_first'])
                         .agg({'allele': lambda x: ','.join(x)})
                         .reset_index())
        # sort to alphabetize
        df_diag_merge.sort_values('allele_first', inplace=True)
        df_diag_merge = df_diag_merge.reset_index()

        df_diag_merge['allele_group_name'] = [x.split('-')[1] if x.count('-') > 1 else x.split('-')[0] for x in
                                              df_diag_merge['allele_first']]
        df_diag_merge['allele_group_name'] = ['_'.join(x.split('_')[:3]) for x in df_diag_merge['allele_group_name']]
        # df_diag_merge['allele_group_name'] = [x.replace('__', '_') for x in df_diag_merge['allele_group_name']]
        df_diag_merge['group_number'] = df_diag_merge.groupby(['allele_group_name'])[
            'allele_group_name'].cumcount().add(1)
        df_diag_merge['header'] = ['{0}-{1}g{2}|{3}'.format(species, x, y, z) if z.count(',') > 0 else j for x, y, z, j
                                   in
                                   zip(df_diag_merge['allele_group_name'],
                                       df_diag_merge['group_number'],
                                       df_diag_merge['allele'],
                                       df_diag_merge['allele_first'])]
    # end rename grouping section
    # add the previous vetted miseq file (this contains novel and links to the haplotype file(s))
    df_diag_2 = fasta_to_df(fasta_path=os.path.join(ref_dir, '{0}.fasta'.format(db_ext[i])),
                            header_name='header',
                            sequence_name='SEQUENCE', as_df=True)
    df_diag = pd.concat([df_diag_merge, df_diag_2], ignore_index=True)
    # df_diag['allele'] = [x.split('|')[0] for x in df_diag['allele']]
    df_diag_dedup = df_diag.drop_duplicates(subset=['SEQUENCE'], keep='first', ignore_index=False)

    df_diag_fasta = df_diag_dedup[['header', 'SEQUENCE']]
    df_diag_fasta.sort_values('header', inplace=True)

    with open(os.path.join(ref_dir, '{0}_long.fasta'.format(db_ext[i])), 'w') as out_file:
        for indx, row in df_diag_fasta.iterrows():
            if len(row['SEQUENCE']) > 130:
                write_sequence(name=row['header'],
                               seq=row['SEQUENCE'],
                               out_file=out_file,
                               replace_illegal_char=True,
                               prefix='')

    df_diag_fasta['header'] = ['{0}-{1}'.format(x.split('|')[0], db_ext[i]) if x.count('|') > 0 else x for x in
                               df_diag_fasta['header']]
    df_diag_fasta['header'] = [x if x.endswith('-{}'.format(db_ext[i])) else '{0}-{1}'.format(x, db_ext[i]) for x in
                               df_diag_fasta['header']]
    with open(os.path.join(ref_dir, '{0}_all.fasta'.format(db_ext[i])), 'w') as out_file:
        for indx, row in df_diag_fasta.iterrows():
            if len(row['SEQUENCE']) > 130:
                write_sequence(name=row['header'].split('|')[0],
                               seq=row['SEQUENCE'],
                               out_file=out_file,
                               replace_illegal_char=True,
                               prefix='')
            else:
                print("Seq shorter than min length:", len(row['SEQUENCE']), row['header'].split('|')[0])
    all_paths.append(os.path.join(ref_dir, '{0}_all.fasta'.format(db_ext[i])))
#############################################
# Concatenate the ipd and miseq fasta files #
#############################################
print('Concatenate the ipd and miseq fasta files')
# pathlist = [ipd_all_filepath, miseq_all_filepath]
with open(all_filepath, 'w') as out_file:
    print(all_filepath)
    for path_i in all_paths:
        fasta_sequences = SeqIO.parse(open(path_i), 'fasta')
        for fasta in fasta_sequences:
            name = fasta.description
            if len(name) > 1:
                write_sequence(name=name,
                               seq=str(fasta.seq),
                               out_file=out_file,
                               replace_illegal_char=True,
                               prefix='')

### #######################################################################
# Create Lookup jsons and create missing list sequences for each database #
###########################################################################
print('Create Lookup jsons and create missing list of sequences for each db')
ipd_fasta_dict = fasta_to_df(fasta_path=all_filepath)
db_dicts = {}
# create a starting blank dictionaries
for suffix in db_ext:
    db_dicts[suffix] = {}
# convert the fasta to keyvalue dictionaries for matching subsets
for k, v in ipd_fasta_dict.items():
    for suffix in db_ext:
        if k.endswith(f'-{suffix}'):
            db_dicts[suffix][k] = v
print(db_ext)
# go through th heiarchy possibilities similar to before, minus the subsets
for i in range(1, len(db_ext)):
    print(db_ext[i])
    compare_suffix = db_ext[i]
    suffix_list = db_ext[0:i]
    for suffix in suffix_list:
        print(db_ext[i], suffix)
        lookup_dict = {}
        for ipd_name, ipd_seq in db_dicts[suffix].items():
            match_list = []
            for diag_name, diag_seq in db_dicts[compare_suffix].items():
                if (diag_seq in ipd_seq) or (rev_comp_st(diag_seq) in ipd_seq):
                    match_list.append(diag_name)
            lookup_dict[ipd_name] = match_list

            # out put to gson
            # ipd_to_exon_lookup_path
            # missing_ipd_alleles_path

        with open(os.path.join(ref_dir, f'{suffix}_to_{compare_suffix}_lookup.json'), 'w') as convert_file:
            convert_file.write(json.dumps(lookup_dict))
        # get the missing list
        missing_ipd_list = []
        for ipd_seq, diag_seq_list in lookup_dict.items():
            if len(diag_seq_list) < 1:
                missing_ipd_list.append(ipd_seq)
        missing_ipd_list.sort()
        # output do data frame
        df_missing_ipd = pd.DataFrame({'allele': missing_ipd_list})
        df_missing_ipd.to_csv(os.path.join(ref_dir, f'{suffix}_to_{compare_suffix}_missing.csv'),
                              header=False,
                              index=False)

        reverse_lookup = {}
        for diag_name, diag_seq in db_dicts[compare_suffix].items():
            match_list = []
            for ipd_name, ipd_seq in db_dicts[suffix].items():
                if (diag_seq in ipd_seq) or (rev_comp_st(diag_seq) in ipd_seq):
                    match_list.append(ipd_name)
            reverse_lookup[diag_name] = match_list

        # out put to gson
        with open(os.path.join(ref_dir, f'{compare_suffix}_to_{suffix}_lookup.json'), 'w') as convert_file:
            convert_file.write(json.dumps(reverse_lookup))
        # get the missing list
        missing_diag_list = []
        for diag_seq, ipd_seq_list in reverse_lookup.items():
            if len(ipd_seq_list) < 1:
                missing_diag_list.append(diag_seq)
        missing_diag_list.sort()
        # output do data frame
        df_missing_diag = pd.DataFrame({'allele': missing_diag_list})
        df_missing_diag.to_csv(os.path.join(ref_dir, f'{compare_suffix}_to_{suffix}_missing.csv'), header=False,
                               index=False)


Path(os.path.join(ref_dir, 'mask_range.csv')).touch()
print('Reference-Files Completed')
if (haplo_fasta is None) or (not os.path.exists(haplo_fasta)):
    print("The haplo type summary files are not specified or do not exist.  They will not be created.")
    exit(0)
###############################
#  Create haplotype Library   #
# - Based on the legacy file  #
###############################
########################################
# open the lookup .json for haplotypes #
# ########################################
print('Create haplotype Library')
with open(haplotype_json_path) as f_in:
    haplo_dict = json.load(f_in)
df = pd.DataFrame(haplo_dict)
df = df.reset_index()
df.rename(columns={'index': 'HAPLOTYPE_CALL'}, inplace=True)
haplotype_cols = list(df.columns)
haplotype_cols = [x for x in haplotype_cols if x.startswith('MHC')]
# stack the data frame so all the columns line up better
df_melt = pd.melt(df, id_vars=['HAPLOTYPE_CALL', 'PREFIX'], value_vars=haplotype_cols,
                  var_name='TYPE', value_name='allele')
# drop the na entries that are created from the dict to df conversion
df_melt.dropna(inplace=True)
df_melt['CALL_COUNT'] = [len(x) for x in df_melt['allele']]
# explode the lists so each list item gets a row
df_melt = df_melt.explode('allele')
# Convert the miseq fasta to a data frame
df_diag_haplo_ref = fasta_to_df(haplo_fasta, header_name='allele', sequence_name='SEQUENCE', as_df=True)
# merge with the data frame by name with the haplotype calling dictionary
df_melt['haplo_allele'] = [like_join_haplo(x, df_diag_haplo_ref, 'allele', 'allele') for x in df_melt['allele']]
df_haplo_call = df_melt.explode('haplo_allele')
df_haplo_call = df_haplo_call[~df_haplo_call['haplo_allele'].isnull()]
df_diag_haplo_ref_m = df_diag_haplo_ref.rename(columns={'allele': 'haplo_allele'})
df_haplo_call_m = df_haplo_call.merge(df_diag_haplo_ref_m, on='haplo_allele')
df_bait_ref = fasta_to_df(fasta_path=all_filepath, as_df=True)
# df_haplo_call_m
# take the reverse complement (rc) of each row and like join to both the forward and rc sequence
print(df_bait_ref)
df_haplo_call_m['bait_allele'] = [like_join_haplo(x, df_bait_ref, 'SEQUENCE', 'allele') for x in
                                  df_haplo_call_m['SEQUENCE']]


df_haplo_call_e = df_haplo_call_m.explode('bait_allele')
df_haplo_call_e.drop(columns=['haplo_allele', 'SEQUENCE'], inplace=True)
df_haplo_call_e.to_csv(haplotype_csv_path, index=False)
