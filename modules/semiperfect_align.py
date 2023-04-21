#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentTypeError
from subprocess import Popen, PIPE
import os


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
    parser.add_argument('--fastq_dir',
                        type=str,
                        help='directory to your input fastq files',
                        required=True)
    parser.add_argument('--bam_dir',
                        type=str,
                        help='directory your output bam files will reside',
                        required=True)
    parser.add_argument('--cp_dir',
                        type=str,
                        help='directory your bbmap executables are',
                        required=True)

    parser.add_argument('--ref_dir',
                        type=str,
                        help='where the config files are stored',
                        required=True)

    parser.add_argument('--bait_fasta',
                        type=str,
                        help='Where you ipd-diag combined fasta exists',
                        default=None,
                        required=False)
    parser.add_argument('--threads',
                        type=int,
                        help='Where you ipd-diag combined fasta exists',
                        default=1,
                        required=False)
    parser.add_argument('--ram',
                        type=int,
                        help='Where you ipd-diag combined fasta exists',
                        default=8000,
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
cp_dir = args.cp_dir
fastq_dir = args.fastq_dir
bam_dir = args.bam_dir
ref_dir = args.ref_dir
bait_fasta = args.bait_fasta
threads = args.threads
ram = args.ram

if bait_fasta is None:
    bait_fasta = os.path.join(ref_dir, 'bait.fasta')

fastq_filelist = os.listdir(fastq_dir)
fastq_filelist = [os.path.join(fastq_dir, x) for x in fastq_filelist if
                  x.endswith('R1_001.fastq.gz') and not x.startswith('._')]
fastq2_filelist = [x.replace('_R1_', '_R2_') for x in fastq_filelist]
sample_list = [os.path.basename(x).split('_')[0] for x in fastq_filelist]
os.makedirs(bam_dir, exist_ok=True)
def print_time(line, update_reads, total_time, samples_procesed):
    try:
        total_time += float(line.strip())
    except:
        pass
    print('{0} seconds to process last {1}M reads. {2:.1f} total time / {3:.1f}M total reads'.format(line.strip(),
                                                                                                     int(update_reads/1000000),
                                                                                                     total_time,
                                                                                                     int(samples_procesed/1000000)))
for in_1, in_2, sample_i in zip(fastq_filelist, fastq2_filelist, sample_list):
    if os.path.exists(os.path.join(bam_dir, '{0}.bam'.format(sample_i))):
        print(os.path.join(bam_dir, '{0}.bam'.format(sample_i)))
        continue
    if not os.path.exists(in_1) and not os.path.exists(in_2):
        print(in_1)
        continue
    update_reads=1000000
    total_time = 0
    cmd_list = ['java' ,
                '-ea',
                '-Xmx{0}m'.format(int(ram)),
                '-Xms{0}m'.format(int(ram)),
                '-cp', '{0}'.format(cp_dir),
                'align2.BBMap', 'build=1',
                'in={0}'.format(in_1),
                'in2={0}'.format(in_2),
                'ref={0}'.format(bait_fasta),
                'outm={0}/{1}.bam'.format(bam_dir,sample_i),
                'semiperfectmode=t',
                'threads={0}'.format(int(threads)),
                'nodisk=f',
                'showprogress2={0}'.format(int(update_reads))]


    print(' '.join(cmd_list))

    run_started = False
    sample_finished = False
    samples_procesed= 0
    with Popen(cmd_list, stdout=PIPE, stderr=PIPE, bufsize=1, universal_newlines=True) as p1:
        for line in iter(p1.stderr.readline, b''):
            #print(line)
            if len(line) == 0:
                break
            check_error = True
            if not run_started:
                print(line.strip())
                if 'mapping threads.' in line.strip():
                    run_started = True
            else:
                samples_procesed += update_reads
                if 'Detecting finished threads' in line.strip():
                    sample_finished = True
                if sample_finished:
                    print(line.strip())
                else:
                    try:
                        total_time += float(line.strip())
                    except:
                        pass
                    print_time(line, update_reads, total_time, samples_procesed)

            checksum_passed = False
        for line in iter(p1.stdout.readline, b''):
            if len(line) == 0:
                break
            print(line.strip())
            checksum_passed = False
