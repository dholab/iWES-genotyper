#!/usr/bin/env python3

"""
TODO
"""

import argparse
import asyncio
import json
import os
from argparse import ArgumentParser, ArgumentTypeError
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@dataclass
class DatabaseRecords:
    """
    TODO
    """

    label: str
    records: Set[SeqRecord]


def parse_command_line_args() -> argparse.Namespace:
    """
    Parses the python script arguments from bash and makes sure files/inputs are valid
    """

    parser = ArgumentParser()
    parser.add_argument(
        "--ref_dir", type=Path, help="where the config files are stored", required=True
    )
    parser.add_argument(
        "-d",
        "--db_path",
        nargs="+",
        help="Database file paths space separated to add multiple -d /path/one /path/two /path/three",
        required=True,
    )
    parser.add_argument(
        "-e",
        "--db_ext",
        nargs="+",
        default=None,
        help="extensions to add at end of the headers to descriminate files (i.e. -e gen exon miseq",
        required=False,
    )
    parser.add_argument(
        "--haplo_fasta",
        type=str,
        help="the database fasta file the haplotype_json is based on",
        default=None,
        required=False,
    )
    parser.add_argument("--species", type=str, help="Species prefix", required=True)
    parser.add_argument(
        "--haplotype_json_path",
        type=str,
        help="haplotype_json_path",
        default=None,
        required=False,
    )
    parser.add_argument(
        "--cp_path",
        type=str,
        help="Directory Where your bbmap executables are stored",
        required=True,
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of threads to run bbmap",
        default=1,
        required=False,
    )
    parser.add_argument(
        "--ram",
        type=int,
        help="Directory ram to dedicate to bbmap java",
        default=8000,
        required=False,
    )
    parser.add_argument(
        "--minimap2_path",
        type=str,
        help="location of minimap2 executable often it is in ./bin and can be called with minimap2",
        default="minimap2",
        required=False,
    )

    return parser.parse_args()


def write_renamed_seq(
    name: str, seq: str, out_file, replace_illegal_char=True, prefix="", suffix=""
):
    """
    TODO
    """
    if replace_illegal_char:
        name = name.replace(
            "*", "_"
        )  # add a double because i want to fine the 02 at the beginning easir
        name = name.replace(":", "_")
    out_file.write(
        f">{prefix}{name}{suffix}\n"
    )  # . We do mostly non human primates I want to differentiate.
    out_file.write(f"{seq}\n")


def fasta_to_df(
    fasta_path=None, header_name="allele", sequence_name="SEQUENCE", as_df=False
):
    """
    TODO
    """
    fasta_sequences = SeqIO.parse(open(fasta_path), "fasta")
    fasta_dict = {}
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        fasta_dict[name] = sequence
    fasta_sequences.close()
    if as_df:
        return pd.DataFrame(fasta_dict.items(), columns=[header_name, sequence_name])
    return fasta_dict


def rename_ref_entries(
    db_exts: List[str], db_paths: List[str], ref_dir: Path
) -> List[DatabaseRecords]:
    """
    Rename and shorten FASTA sequence identifiers in reference databases.

    Args:
        db_exts (List[str]): List of database extension labels.
        db_paths (List[str]): List of paths to FASTA database files.
        ref_dir (Path): Directory containing input FASTA files.

    Returns:
        Dict[str, SeqRecord.SeqRecord]: Dictionary mapping database labels to
        sets of renamed SeqRecord objects.

    Raises:
        AssertionError: If number of db_exts and db_paths don't match,
                        if ref_dir doesn't exist, or if the number of
                        successfully renamed databases doesn't match inputs.
    """
    # make sure the number of databases and extensions match
    assert len(db_paths) == len(db_exts), f"""
    Please double check that the provided reference database locations:
    
    {db_paths}
    
    ...equals the number of database extensions:
    
    {db_exts}
    """

    # make sure the reference directory exists
    assert os.path.isdir(ref_dir), f"""
    Please double check that the provided reference directory {ref_dir} exists.
    """
    # TODO - switch to proper logging
    print(
        f"""
    Databases mapped to their labels in this way:

    {zip(db_exts, db_paths)}
    """
    )

    renamed_dbs: List[DatabaseRecords] = []
    for label, path in zip(db_exts, db_paths):
        # pull in a generator for the FASTA
        ref_sequences = SeqIO.parse(path, "fasta")

        # use a set comprehension to create a new set of SeqRecord objects with names
        # abbreviated to the first text from the left up until the first pipe '|'
        renamed_records = {
            SeqRecord(record.seq, id=record.name.split("|")[0])
            for record in ref_sequences
        }
        bundled_records = DatabaseRecords(label, records=renamed_records)

        renamed_dbs.append(bundled_records)

    assert len(renamed_dbs) == len(db_exts), f"""
    The length of the successfully renamed databases, {len(renamed_dbs)}, does not match the
    length of the provided input databases, {len(db_exts)}.
    """

    # TODO - switch to proper logging
    print(f"Successfully shortened names for {len(renamed_dbs)} reference FASTAs.")

    return renamed_dbs


def find_ref_overlaps(
    ref_dict: List[DatabaseRecords], concat_seqs: List[SeqRecord] = []
):
    """
    TODO
    """

    # If we've
    if len(ref_dict) < 2:
        return concat_seqs


def main() -> None:
    """
    TODO
    """
    args = parse_command_line_args()

    # make sure that the same number of database paths and extensions were provided
    assert len(args.db_path) == len(args.db_ext) and len(args.db_path) > 0, f"""
    Please double check that the provided reference database locations:
    
    {args.db_path}
    
    ...equals the number of database extensions:
    
    {args.db_ext}
    """

    # make sure the reference directory exists
    assert os.path.isdir(args.ref_dir), f"""
    Please double check that the provided reference directory {args.ref_dir} exists.
    """

    # compute a new set of renamed records for each input database to maintain
    # compatability for some the tools used downstream
    renamed_ref_dict: Dict[str, SeqRecord] = rename_ref_entries(
        args.db_ext, args.db_paths, args.refdir
    )

    haplotype_json_path = args.haplotype_json_path
    haplo_fasta = args.haplo_fasta
    cp_path = args.cp_path
    threads = args.threads
    ram = args.ram

    species = args.species
    minimap2_path = args.minimap2_path

    # TODO - replace with proper logging
    print(minimap2_path)

    split_char_ipd = " "
    split_char_diag = ""

    os.makedirs(args.ref_dir, exist_ok=True)

    all_filepath = os.path.join(args.ref_dir, "bait.fasta")
    haplotype_csv_path = os.path.join(args.ref_dir, "haplotype_lookup.csv")
    pd.options.mode.chained_assignment = None


if __name__ == "__main__":
    main()
