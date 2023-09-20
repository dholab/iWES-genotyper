import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from os.path import join as join_path
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
    parser.add_argument('--in_path',
                        type=str,
                        help='file path to the distribution csv generated form mapped_read_size_dist.py',
                        required=True)
    parser.add_argument('--out_dir',
                        type=str,
                        help='path to your output dir.',
                        required=True)
    parser.add_argument('--prefix',
                        type=str,
                        help='path to your output dir.',
                        required=True)
    parser.add_argument('--sample',
                        type=int,
                        help='path to your output dir.',
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
in_path = args.in_path
prefix = args.prefix
sample = args.sample
# Apply the default theme
# sns.set_theme()
print(f'Opening file {in_path}')
df_plot = pd.read_csv(in_path)
df_plot_i = df_plot[df_plot['gs_id'] == sample]
df_plot_i = df_plot_i[df_plot_i['ALLELE'] != 'Mamu-A1_028_01_01_01-gen']
# df_plot_i.to_csv(join_path(out_dir, f'{prefix}_{sample}_depth_of_coverage_plot.csv'), index=False)
df_plot_i = df_plot_i[['COMPARING_POSITION', 'DEPTH_ADJ', 'DEPTH', 'ALLELE']]
# df_plot_i.rename(columns={'COMPARING_POSITION': 'POSITION', 'DEPTH_ADJ': 'DEPTH'})

sns.lineplot(data=df_plot_i, x="COMPARING_POSITION", y="DEPTH", hue="ALLELE")
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.ylim((-1,185))
plt.gcf().set_size_inches(10, 5)
print('Saving plot')
savepath=join_path(out_dir, f'{prefix}_{sample}_depth_of_coverage_plot_noadj.svg')
plt.savefig(savepath, dpi=200)
print(f'Completed plot: {savepath}')
