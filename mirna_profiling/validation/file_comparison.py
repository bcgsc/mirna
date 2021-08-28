#!/usr/bin/env python

import argparse
import os
import math
import yaml
import numpy
from pathlib import Path
import json
import pandas
import gzip

# from filecmp import dircmp

REL_TOL = 1e-4


def main():
    """
    A script to compare all files in two directories.
    Used to validate outputs of a new pipeline against those of an existing pipeline.
    """
    args = parse_args()
    if args.file_types:
        with open(args.file_types, "r") as f:
            file_types = yaml.safe_load(f)
    else:
        file_types = None

    compare_directory(args.dir1, args.dir2, file_types=file_types)


def compare_directory(dir1, dir2, file_types=None):
    """

    Args:
        dir1 (str):
        dir2 (str):

    Returns:

    """
    files = Path(dir1).rglob("*")
    are_identical = True
    for infile in files:
        infile = infile.resolve()
        if infile.is_dir():
            continue

        outfile = str(infile).replace(os.path.abspath(dir1), os.path.abspath(dir2))
        if not os.path.exists(outfile):
            verdict = "missing"
            are_identical = False
        elif not compare_files(infile, outfile, file_types=file_types):
            are_identical = False
            verdict = "different"
        else:
            verdict = "identical"
        print("{} : {}".format(infile, verdict))

    return are_identical


def compare_files(file_1, file_2, file_types=None, **kwargs):

    fname = os.path.basename(file_1)
    if fname.endswith(".csv") or (file_types and fname in file_types.get("csv")):
        are_identical = compare_csv(file_1, file_2, **kwargs)
    elif fname.endswith(".json"):
        are_identical = compare_json(file_1, file_2, **kwargs)
    elif fname.endswith(".gz"):
        are_identical = compare_gzip(file_1, file_2, **kwargs)
    else:
        are_identical = compare_text(file_1, file_2, **kwargs)

    return are_identical


def compare_csv(file_1, file_2, sep=None, index_col=None):
    are_identical = False
    with open(file_1, "r") as f:
        first_line = f.readline().strip()
    if not sep:
        sep = "\t" if "\t" in first_line else ","
    if index_col is not None:
        df1 = pandas.read_csv(
            file_1, index_col=index_col, header=0, sep=sep
        ).sort_index()
        df2 = pandas.read_csv(
            file_2, index_col=index_col, header=0, sep=sep
        ).sort_index()
    else:
        df1 = pandas.read_csv(file_1, header=0, sep=sep)
        df2 = pandas.read_csv(file_2, header=0, sep=sep)

    if set(df1.columns) == set(df2.columns):
        txt_cols = [col for col in df1.columns if df1[col].dtype == "object"]
        num_cols = [col for col in df1.columns if df1[col].dtype != "object"]

        if all(
            set(numpy.isclose(df1[x], df2[x], rtol=REL_TOL).tolist()) == {True}
            for x in num_cols
        ) and df1[txt_cols].equals(df2[txt_cols]):
            are_identical = True

    return are_identical


def compare_gzip(file_1, file_2, ordered=False):
    with gzip.open(file_1, "rb") as file1:
        with gzip.open(file_2, "rb") as file2:
            if not ordered:
                are_identical = sorted(file1) == sorted(file2)
            else:
                are_identical = file1 == file2

    return are_identical


def compare_text(file_1, file_2, ordered=False):
    with open(file_1, "r") as file1:
        with open(file_2, "r") as file2:
            if not ordered:
                are_identical = sorted(file1) == sorted(file2)
            else:
                are_identical = file1 == file2

    return are_identical


def compare_json(file_1, file_2, ordered=False, key=None):
    with open(file_1, "r") as file1:
        with open(file_2, "r") as file2:
            records1 = json.load(file1)
            records2 = json.load(file2)
    if not ordered:
        records1 = sorted(records1, key=key)
        records2 = sorted(records2, key=key)
    are_identical = len(records1) == len(records2) and all(
        x == y for x, y in zip(records1, records2)
    )

    return are_identical


def parse_args():
    parser = argparse.ArgumentParser(description="Compare all files in two directories")
    parser.add_argument("--dir1", help="One directory path", required=True)
    parser.add_argument("--dir2", help="The other directory path", required=True)
    parser.add_argument(
        "--file_types",
        help="file_types yaml",
        required=True,
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
