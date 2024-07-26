#!/usr/bin/env python3

"""
TODO
"""

import json
import os
from argparse import ArgumentTypeError

from Bio import SeqIO


def str2bool(v):
    """
    TODO
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise ArgumentTypeError("Boolean value expected.")


def generate_revcomp(seq: str) -> str:
    """
    TODO
    """
    assert (
        len(seq) > 0
    ), f"The provided sequence, {seq}, appears to be empty or otherwise corrupted."
    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    seq = seq[::-1]
    return seq


def trim_read_for_diag(seq: str, cigartupple):
    """
    TODO
    """
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


def single_file_per_fasta(fasta_path, single_fasta_dir, name_num_json_path):
    i = 0
    ipd_num = {}
    fasta_sequences = SeqIO.parse(open(fasta_path), "fasta")
    for fasta in fasta_sequences:
        i += 1
        name, sequence = fasta.id, str(fasta.seq)
        if len(sequence) > 1:
            truncated_ref = os.path.join(single_fasta_dir, "{0}.fasta".format(i))
            ipd_num[name] = "{0}".format(i)
            with open(truncated_ref, "w") as out_file:
                out_file.write(">{0}\n".format(name))
                out_file.write("{0}\n".format(sequence))
    with open(name_num_json_path, "w") as convert_file:
        convert_file.write(json.dumps(ipd_num))
    return


def like_join(x, df_ref, column_i="SEQUENCE"):
    """
    TODO
    """
    name_list = []
    for idx, row in df_ref.iterrows():
        if x in row["allele"]:
            name_list.append(row[column_i])
    return name_list


def rev_comp_like_join(inner_seq, df_ref):
    name_list = []
    for idx, row in df_ref.iterrows():
        if (inner_seq in row["SEQUENCE"]) or (
            generate_revcomp(inner_seq) in row["SEQUENCE"]
        ):
            name_list.append(row["allele"])
    return name_list
