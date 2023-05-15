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

# Apply the default theme
# sns.set_theme()
print(f'Opening file {in_path}')
df_plot = pd.read_csv(in_path)

print('Generating plot')
sns.displot(df_plot, x="MAPPED_LENGTH", kind="kde", bw_adjust=2)
plt.gcf().set_size_inches(5, 5)
print('Saving plot')
savepath=join_path(out_dir, 'read_length_dist.svg')
plt.savefig(savepath, dpi=200)
print(f'Completed plot: {savepath}')