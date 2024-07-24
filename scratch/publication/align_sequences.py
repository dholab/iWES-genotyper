import pandas as pd
import os
from Bio import SeqIO
import os
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
from argparse import ArgumentParser, ArgumentTypeError
import itertools

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
                        help='filepath to the fasta. The fasta seq_ids must match ids in the  protein groups path',
                        required=True)
    parser.add_argument('--depth_dir',
                        type=str,
                        help='path to ypur depth files.',
                        required=True)
    parser.add_argument('--out_dir',
                        type=str,
                        default='./out',
                        help='Directory files are outputted',
                        required=False)
    parser.add_argument('--peptide_length',
                        type=int,
                        default=1,
                        help='Set to the maximum Peptide length in your seq_path',
                        required=False)
    parser.add_argument('--comparing_strain_prefix',
                        type=str,
                        default="COMPARE",
                        help='Uses the SEQ_NAME (first column) of prot_groups_path table if blank, prefiex assigned in column names for csv files',
                        required=False)
    parser.add_argument('--muscle_path',
                        type=str,
                        default="/muscle",
                        help='path to muscle program, https://www.drive5.com/muscle',
                        required=False)
    parser.add_argument('--docker_image',
                        type=str,
                        default="pepmeld:v1_1",
                        help='path to docker image',
                        required=False)
    parser.add_argument('--use_docker',
                        type=str,
                        default=False,
                        help='path to docker image',
                        required=False)
    parser.add_argument('--overwrite',
                        type=str2bool,
                        default=False,
                        required=False
                        )
    parser.add_argument('-c',
                        '--mhc_class',
                        nargs='+',
                        default=['Mamu-A1', 'Mamu-B'],
                        help='extensions to add at end of the headers to descriminate files (i.e. -e gen exon miseq',
                        required=False)
    parser.add_argument('-e',
                        '--db_ext',
                        nargs='+',
                        default=['gen'],
                        help='extensions to add at end of the headers to descriminate files (i.e. -e gen exon miseq',
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
# same arguments to a local variable by same name as the argument
ref_dir = args.ref_dir
out_dir = args.out_dir
peptide_length = args.peptide_length
muscle_path = args.muscle_path
docker_image = args.docker_image
overwrite = args.overwrite
mhc_class_list = args.mhc_class
db_ext_list = args.db_ext
depth_dir = args.depth_dir
use_docker = args.use_docker

fasta_path = os.path.join(ref_dir, 'bait.fasta')
ref_out_dir = os.path.join(ref_dir, 'aligned')
os.makedirs(out_dir, exist_ok=True)
os.makedirs(ref_out_dir, exist_ok=True)


def get_aligned_positions_df(df_fasta,
                             seq_comp_id,
                             comparing_strain_prefix="COMPARING",
                             peptide_length=1):
    import pandas as pd
    import math
    # Easier to work with dictionary than iterrate the rows of data frame
    seq_dict = df_fasta.set_index('SEQ_ID')['SEQUENCE'].to_dict()
    seq_comp = seq_dict[seq_comp_id]
    seq_aligned_positions = {}
    total_aligned_score = {}
    seq_alignment_score_list = {}
    df_all = pd.DataFrame()
    for seq_id, seq in seq_dict.items():
        if seq == '':
            continue
        alignment_score = []
        # First over arching position.  Tracks the position of the comparing sequence
        j = 0
        score = 0
        # uses this to track the comparing sequence position
        comp_position = []
        # loop through the length of the sequence
        # Because the file (should) be aligned they should be the same length
        for i in range(0, len(seq)):
            # add to j because it is not a -
            if seq_comp[i] != '-':
                j = j + 1

            if seq[i] != '-':
                # this is the position that is matched to the comparing seq
                comp_position.append(j)
                # add to the total socre if they are the same (already check for -)
                if seq[i] == seq_comp[i]:
                    score = score + 1
                # track the end of the comparing sequence as it loops to i + k
                k = i
                # used to track valid positions of the comparing sequence
                m = 0
                # used to track the valid positions of the sequece of interest in this loop
                p = 0
                align_score = 0
                # Seq length
                # loops to the end of the peptide.
                # p and m are tracked for non - (amino acids declared)
                while (p < peptide_length) and (m < peptide_length):
                    # Break if the comparing sequence is at the end.
                    # Nothing else to align
                    if k + 1 > len(seq_comp):
                        break
                    # Track the comparing sequence non -  value count
                    # vs peptide length
                    if seq_comp[k] != '-':
                        m = m + 1
                    # Track the sequence non -  value count
                    # vs peptide length
                    if seq[k] != '-':
                        p = p + 1
                    # If both are not - and match add one to the align score
                    if (seq_comp[k] != '-') and (seq[k] != '-'):
                        if seq_comp[k] == seq[k]:
                            align_score = align_score + 1

                    k = k + 1

                alignment_score.append(align_score)

        # add to the dictionaries to store the values by seq_Id
        # seq_aligned_positions is the position of the comparing sequence
        seq_aligned_positions[seq_id] = comp_position
        # total_aligned_score is the total match count (subtracting '-' for both)
        # useful in ordering the sequences by over all mismatch count
        total_aligned_score[seq_id] = score
        # this is a list of each alignment match count
        seq_alignment_score_list[seq_id] = alignment_score
        df = pd.DataFrame({'SEQ_ID': seq_id,
                           'POSITION': list(range(1, len(comp_position) + 1)),
                           '{0}_POSITION'.format(comparing_strain_prefix): comp_position,
                           'MATCH_COUNT': alignment_score,
                           'SCORE': score})
        df_all = pd.concat([df_all, df], ignore_index=True)
    return df_all


def open_fasta_as_df(fasta_path):
    with open(fasta_path) as fasta_file:  # Will close handle cleanly
        seq_id = []
        mhc_class = []
        seq = []
        db = []
        for title, sequence in SimpleFastaParser(fasta_file):
            seq_id.append(title)  # First word is ID
            mhc_class.append(title.split('_')[0])
            seq.append(sequence)
            db.append(title.split('-')[-1])
    return pd.DataFrame({'SEQ_ID': seq_id, 'MHC_CLASS': mhc_class, 'SEQUENCE': seq, 'DB': db})


# same arguments to a local variable by same name as the argument
depth_all_path = os.path.join(out_dir, 'mhc_class_depth.csv.gz')
if not os.path.exists(depth_all_path) or overwrite:
    file_list = os.listdir(depth_dir)
    file_list = [os.path.join(depth_dir,x) for x in file_list if not x.startswith('._') and x.endswith('depth-gen.csv')]
    df_depth = pd.DataFrame()
    i = 0
    for filepath_i in file_list:
        i +=1
        print(filepath_i, i, len(file_list))
        df_depth_i = pd.read_csv(filepath_i)
        df_depth_i.rename(columns={'SAMPLE_NUM':'gs_id'}, inplace=True)
        df_depth_i = df_depth_i[df_depth_i['HAPLOTYPE_GROUP'].isin(mhc_class_list)]
        df_depth = pd.concat([df_depth, df_depth_i], ignore_index=True)
    print(f'Saving {depth_all_path}')
    df_depth.to_csv(depth_all_path, index=False)
else:
    print(f'Opening {depth_all_path}')
    df_depth = pd.read_csv(depth_all_path)
print(f'Opening {fasta_path}')
df_fasta = open_fasta_as_df(fasta_path)

for mhc_class_i, db_i in itertools.product(mhc_class_list, db_ext_list):
    print(mhc_class_i, db_i)
    df_fasta_i = df_fasta[(df_fasta['MHC_CLASS'] == mhc_class_i) & (df_fasta['DB'] == db_i)]
    header_list = list(df_fasta_i['SEQ_ID'].unique())
    out_path = os.path.join(ref_out_dir, '{0}__{1}_aligned.fasta'.format(mhc_class_i, db_i))
    in_path = os.path.join(ref_out_dir, '{0}__{1}_filtered.fasta'.format(mhc_class_i, db_i))
    # create a new fasta with only from teh columns
    print(in_path)
    with open(in_path, "w") as f:
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.description in header_list:
                f.write(record.format("fasta"))
    if not os.path.exists(out_path) or overwrite:
        print(f'Running muscle on {in_path} to {out_path}')
        # docker run -v /Volumes/T8:/Volumes/T8 pepmeld:v1_1 /muscle -in {in_path} -out {out_path}
        if use_docker:
            print('using docker')
            subprocess.call(['docker',
                             'run',
                             '-v',
                             f'{in_path}:/in',
                             '-v',
                             f'{out_path}:/out',
                             docker_image,
                             muscle_path,
                             '-in',
                             '/in',
                             '-out',
                             '/out'],
                            shell=False)
        else:
            print('using local muscle')
            subprocess.call([muscle_path,
                             '-in',
                             in_path,
                             '-out',
                             out_path],
                            shell=False)
    aligned_csv_path = os.path.join(out_dir, f'{mhc_class_i}__{db_i}_aligned_fasta.csv.gz')
    if not os.path.exists(aligned_csv_path) or overwrite:
        print(f'Aggregating and Saving aligned csv to {aligned_csv_path}')
        df_aligned_fasta = open_fasta_as_df(out_path)
        df_aligned_fasta.sort_values(by='SEQ_ID', inplace=True)
        seq_comp_id = list(df_aligned_fasta['SEQ_ID'])[0]

        df_aligned = get_aligned_positions_df(df_fasta=df_aligned_fasta,
                                              seq_comp_id=seq_comp_id,
                                              comparing_strain_prefix='COMPARING',
                                              peptide_length=peptide_length)

        df_aligned.to_csv(aligned_csv_path, index=False)
    else:
        print(f'Opening aligned csv: {aligned_csv_path}')
        df_aligned = pd.read_csv(aligned_csv_path)


    print(f'Merging and aggregating aligned csv to aggregated depth files')
    df_aligned.rename(columns={'SEQ_ID': 'ALLELE'}, inplace=True, errors='ignore')
    df_aligned_m = df_aligned.merge(df_depth,
                                    on=['POSITION', 'ALLELE'],
                                    how='inner')
    aligned_depth_path = os.path.join(out_dir, f'aligned_depth_{mhc_class_i}__{db_i}.csv.gz')
    print(f'Saving aligened depth file, not aggregated by compared position: {aligned_depth_path}')
    df_aligned_m.to_csv(aligned_depth_path, index=False)
    agg_mean = ['POSITION',
                'unique_maps_per_allele',
                'DEPTH',
                'num_read_maps',
                'DEPTH_ADJ',
                'DEPTH_RATIO']
    df_aligned_m_agg = df_aligned_m.groupby(['COMPARING_POSITION', 'ALLELE', 'gs_id'], as_index=False)[agg_mean].mean()
    df_aligned_m_agg[agg_mean] = df_aligned_m_agg[agg_mean].round(2)
    aligned_depth_agg_path = os.path.join(out_dir, f'aligned_depth_{mhc_class_i}__{db_i}_aggregated.csv')
    print(f'Saving aligned depth file: aggregating aligned csv to aggregated depth files to: {aligned_depth_agg_path}')
    df_aligned_m_agg.to_csv(aligned_depth_agg_path, index=False)
