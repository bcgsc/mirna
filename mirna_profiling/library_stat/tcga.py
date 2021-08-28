#!/usr/bin/env python

import argparse
import logging
import os
from pathlib import Path

import pandas

from mirna_profiling.annotation.mirbase import Mirbase

NORM_FACTOR = 1000000


def main():
    args = parse_args()

    mirbase = Mirbase(args.mirbase, args.species_code)
    for sample_dir in (str(x.parent) for x in Path(args.proj_dir).rglob("*_features")):
        write_tcga_files(sample_dir, args.ucsc, mirbase, out_dir=args.out_dir)


def write_tcga_files(sample_dir, genome, mirbase, out_dir=None):
    """
    Generate TCGA-formatted expression quantification mirnas.txt and isoforms.txt files
    Depends on miRNA.txt, crossmapped.txt and isoforms.txt output files from alignment_stat

    Args:
        mirbase (Mirbase): an instance of Mirbase class
        genome (str): genome name, e.g. hg38
        sample_dir (str):sample directory path with miRNA expression files under <sample>_features/
        out_dir (str): output file directory

    Returns:
        (bool): True if write is successful, False otherise

    """
    feature_dir = [str(x) for x in Path(sample_dir).glob("*_features")]
    if not feature_dir:
        logging.warning("Feature directory is missing under {}".format(sample_dir))
        return None, None
    feature_dir = feature_dir[0]
    mirna_file = os.path.join(feature_dir, "miRNA.txt")
    isoform_file = os.path.join(feature_dir, "isoforms.txt")
    for mfile in (mirna_file, isoform_file):
        if not os.path.exists(mfile) or os.path.getsize(mfile) == 0:
            logging.warning(
                "File {} is empty. Cannot generate TCGA-formatted file".format(mfile)
            )
            return None, None

    mirna_ids = mirbase.get_mirna_ids()
    gene_df = pandas.DataFrame({"miRNA_ID": mirna_ids})
    gene_df.set_index("miRNA_ID", inplace=True)

    sample_name = os.path.basename(feature_dir).replace("_features", "")
    if not out_dir:
        tcga_dir = os.path.join(feature_dir, "tcga")
    else:
        tcga_dir = os.path.join(out_dir, sample_name)
    os.makedirs(tcga_dir, 0o775, exist_ok=True)
    tcga_mirna = os.path.join(tcga_dir, "mirnas.txt")
    tcga_isoform = os.path.join(tcga_dir, "isoforms.txt")

    df_mirna = pandas.read_csv(
        mirna_file,
        names=["miRNA_ID", "read_count"],
        usecols=[0, 1],
        sep=" ",
        dtype={"read_count": "int64"},
    )
    df_mirna["cross-mapped"] = False

    crossmapped_file = os.path.join(feature_dir, "crossmapped.txt")
    if os.path.exists(crossmapped_file) and os.path.getsize(crossmapped_file) > 0:
        df_cross = pandas.read_csv(
            crossmapped_file,
            names=["miRNA_ID", "read_count"],
            usecols=[0, 1],
            sep=" ",
            dtype={"read_count": "int64"},
        )
        df_cross["miRNA_ID"] = df_cross["miRNA_ID"].apply(lambda x: x.split(";"))
        df_cross["cross-mapped"] = True
        df_cross = df_cross.explode("miRNA_ID", ignore_index=True)
        df_mirna = pandas.concat([df_mirna, df_cross], axis=0)
    df_mirna["miRNA_ID"] = df_mirna["miRNA_ID"].apply(lambda x: x.split(",")[0])
    df_mirna = df_mirna.groupby("miRNA_ID").sum()

    df_mirna = pandas.concat([gene_df, df_mirna], axis=1).fillna(0).astype(int)
    total_mirna = df_mirna["read_count"].sum()
    df_mirna["reads_per_million_miRNA_mapped"] = (
        df_mirna["read_count"] * NORM_FACTOR / total_mirna
    )
    df_mirna["cross-mapped"] = df_mirna["cross-mapped"].apply(
        lambda x: "Y" if x > 0 else "N"
    )
    df_mirna.sort_index(inplace=True)
    df_mirna.to_csv(
        tcga_mirna,
        sep="\t",
        float_format="%0.6f",
        columns=["read_count", "reads_per_million_miRNA_mapped", "cross-mapped"],
    )

    df_isoform = pandas.read_csv(
        isoform_file,
        sep=" ",
        names=[
            "miRNA_ID",
            "chr",
            "start",
            "end",
            "strand",
            "seq",
            "read_count",
            "cross-mapped",
            "miRNA_region",
            "isomir_name",
        ],
    )
    df_isoform["isoform_coords"] = df_isoform.apply(
        lambda x: "{}:{}:{}-{}:{}".format(
            genome, x["chr"], x["start"], x["end"], x["strand"]
        ),
        axis=1,
    )
    df_isoform.sort_values(["miRNA_ID", "chr", "start", "end"], inplace=True)
    df_isoform["reads_per_million_miRNA_mapped"] = (
        df_isoform["read_count"] * NORM_FACTOR / total_mirna
    )
    df_isoform["cross-mapped"] = df_isoform["cross-mapped"].apply(
        lambda x: "Y" if x > 0 else "N"
    )
    df_isoform.to_csv(
        tcga_isoform,
        index=False,
        sep="\t",
        float_format="%0.6f",
        columns=[
            "miRNA_ID",
            "isoform_coords",
            "read_count",
            "reads_per_million_miRNA_mapped",
            "cross-mapped",
            "miRNA_region",
        ],
    )
    return df_mirna, df_isoform


def parse_args():
    parser = argparse.ArgumentParser(
        description="miRNA Gene-by-Sample Expression Matrix"
    )
    parser.add_argument(
        "-p",
        "--proj_dir",
        help="project directory path with annotated b/sam files",
        required=True,
    )

    parser.add_argument(
        "-u",
        "--ucsc",
        help="UCSC genome / database name, e.g. hg38",
        required=True,
    )

    parser.add_argument(
        "-m",
        "--mirbase",
        help="miRBase database version, e.g. mirna_21",
        required=True,
    )

    parser.add_argument(
        "-s",
        "--species_code",
        help="miRBase species code, e.g. hsa for human",
        required=True,
    )

    parser.add_argument(
        "-o",
        "--out_dir",
        help="output directory path",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
