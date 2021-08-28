#!/usr/bin/env python

import argparse
import logging
import numpy as np
import pandas
import os
import warnings
from pathlib import Path

from mirna_profiling.annotation.mirbase import Mirbase, get_common_id

NORM_FACTOR = 1000000


def main():
    """
    Creates miRNA gene-by-sample matrix from miRNA output files from alignment_stat

    """

    args = parse_args()
    mirbase = Mirbase(
        args.mirbase_version, args.species_code, config_file=args.config_file
    )
    alignment_stat_to_exp_matrix(
        args.proj_dir,
        mirbase,
        out_dir=args.out_dir,
    )


def alignment_stat_to_exp_matrix(proj_dir, mirbase, out_dir=None):
    """

    Args:
        mirbase (Mirbase): an instance of Mirbase class
        proj_dir (str): project directory path with annotated b/sam files
        out_dir (str): output directory for expression matrix files

    Returns:

    """
    if not out_dir:
        out_dir = proj_dir
    exp_df = get_mirna_expression_matrix(proj_dir, mirbase)
    exp_to_csv(exp_df, "expn_matrix", out_dir)
    mimat_exp_df = get_mimat_expression_matrix(proj_dir, mirbase)
    exp_to_csv(mimat_exp_df, "expn_matrix_mimat", out_dir)

    if len(exp_df) == len(mirbase.get_mirna_ids()) and len(mimat_exp_df) == len(
        mirbase.get_mature_mirna_acc()
    ):
        matrix_status = True
    else:
        matrix_status = False

    status_file = os.path.join(
        out_dir,
        "expression_matrix.{}".format("success" if matrix_status else "failure"),
    )
    Path(status_file).touch()


def get_mirna_expression_matrix(proj_dir, mirbase):
    """
    Get miRNA expression matrix for project
    Depends on mirna_species.txt files in <sample>_features directory under proj_dir

    Args:
        mirbase (Mirbase): an instance of Mirbase class
        proj_dir (str): project directory path with annotated b/sam files

    """
    logging.info("Generating miRNA expression matrix")
    mirna_ids = mirbase.get_mirna_ids()

    exp_df = pandas.DataFrame({"Gene": mirna_ids})
    exp_df.set_index("Gene", inplace=True)

    proj_dir = os.path.abspath(proj_dir)
    for feature_dir in Path(proj_dir).rglob("*_features"):
        expr_file = os.path.join(feature_dir, "mirna_species.txt")
        sample_name = os.path.basename(feature_dir).split("_")[0].split(".")[0]
        if not os.path.exists(expr_file):
            logging.warning(
                "Missing mirna_species.txt for sample {} in {}".format(
                    sample_name, feature_dir
                )
            )
            continue
        if sample_name in exp_df.columns:
            logging.warning(
                "Duplicate sample {} in project directory {}".format(
                    sample_name, proj_dir
                )
            )
            continue

        df = pandas.read_csv(
            expr_file,
            names=["Gene", sample_name],
            usecols=[0, 1],
            sep=" ",
            dtype={sample_name: "int64"},
        )
        df.set_index("Gene", inplace=True)
        exp_df = pandas.concat([exp_df, df], axis=1)

    exp_df = exp_df.fillna(0).astype(int)
    exp_df = exp_df[sorted(exp_df.columns)]

    return exp_df


def get_mimat_expression_matrix(proj_dir, mirbase):
    """
    Get mature miRNA expression matrix for project

    Depends on miRNA.txt and crossmapped.txt files in <sample>_features directory under proj_dir

    Args:
        mirbase (Mirbase): an instance of Mirbase class
        proj_dir (str): project directory path with annotated b/sam files

    """
    logging.info("Generating mature miRNA expression matrix")
    proj_dir = os.path.abspath(proj_dir)
    gene_df = pandas.DataFrame(mirbase.get_mirna_mimat_ids())
    gene_df["Gene"] = gene_df.apply(
        lambda x: "{}.{}".format(
            get_common_id(
                gene_df.loc[gene_df["mature_acc"] == x["mature_acc"]][
                    "mirna_id"
                ].to_list()
            ),
            x["mature_acc"],
        ),
        axis=1,
    )
    gene_df.drop("mirna_id", axis=1, inplace=True)
    gene_df.drop_duplicates("mature_acc", ignore_index=True, inplace=True)
    gene_df.set_index("mature_acc", inplace=True)

    exp_df = None
    for feature_dir in Path(proj_dir).rglob("*_features*"):
        mirna_file = os.path.join(feature_dir, "miRNA.txt")
        sample_name = os.path.basename(feature_dir).split("_")[0].split(".")[0]
        if not os.path.exists(mirna_file):
            logging.warning(
                "Missing miRNA.txt for sample {} in {}".format(sample_name, feature_dir)
            )
            continue

        if exp_df is not None and sample_name in exp_df.columns:
            logging.warning(
                "Duplicate sample {} in project directory {}".format(
                    sample_name, proj_dir
                )
            )
            continue
        df_mirna = pandas.read_csv(
            mirna_file,
            names=["Gene", sample_name],
            usecols=[0, 1],
            sep=" ",
            dtype={sample_name: "int64"},
        )
        crossmapped_file = os.path.join(feature_dir, "crossmapped.txt")
        if os.path.exists(crossmapped_file) and os.path.getsize(crossmapped_file) > 0:
            df_cross = pandas.read_csv(
                crossmapped_file,
                names=["Gene", sample_name],
                usecols=[0, 1],
                sep=" ",
                dtype={sample_name: "int64"},
            )
            df_cross["Gene"] = df_cross["Gene"].apply(lambda x: x.split(";"))
            df_cross = df_cross.explode("Gene", ignore_index=True)

            df_mirna = pandas.concat([df_mirna, df_cross], axis=0)

        df_mirna = df_mirna.loc[
            df_mirna["Gene"].str.contains("mature")
            | df_mirna["Gene"].str.contains("star")
        ]
        df_mirna["Gene"] = df_mirna.apply(
            lambda x: extract_gene_id(gene_df, x),
            axis=1,
        )
        df_mirna = df_mirna.groupby("Gene").sum()

        if exp_df is None:
            exp_df = df_mirna
        else:
            exp_df = pandas.concat([exp_df, df_mirna], axis=1)
    gene_df = gene_df.reset_index(drop=True).set_index("Gene")
    exp_df = pandas.concat([gene_df, exp_df], axis=1).fillna(0).astype(int)
    exp_df = exp_df[sorted(exp_df.columns)]

    return exp_df


def extract_gene_id(gene_df, mirna):
    acc = mirna["Gene"].split(",")[2]
    try:
        gene = gene_df.at[acc, "Gene"]
    except Exception:
        logging.error(
            "mature miRNA {} in expression matrix file is not present in mirBase. Please check that the same miRBase version is used to generate the matrix file".format(
                acc
            )
        )
        gene = "NA"

    return gene


def exp_to_csv(exp_df, prefix, out_dir):
    """
    Write expression matrix to CSV files

    Args:
        exp_df (pandas.DataFrame):
        prefix (str): prefix for expression CSV files
        out_dir (str): output directory for CSV files

    """
    sorted_columns = sorted(exp_df.columns)
    exp_df.sort_index(inplace=True)
    exp_df.to_csv(
        os.path.join(out_dir, "{}.txt".format(prefix)),
        sep="\t",
        header=True,
        columns=sorted_columns,
        index=True,
    )
    exp_norm_df = exp_df.apply(lambda x: x * NORM_FACTOR / exp_df[x.name].sum(), axis=0)
    exp_norm_df.sort_index(inplace=True)
    exp_norm_df.to_csv(
        os.path.join(out_dir, "{}_norm.txt".format(prefix)),
        float_format="%0.6f",
        sep="\t",
        header=True,
        columns=sorted_columns,
        index=True,
    )

    # ignore warning message
    warnings.filterwarnings("ignore", message="divide by zero encountered")
    exp_norm_log_df = exp_norm_df.apply(lambda x: np.log(x) / np.log(2), axis=0)
    exp_norm_log_df.replace([np.inf, -np.inf], 0, inplace=True)
    exp_norm_log_df.sort_index(inplace=True)
    exp_norm_log_df.to_csv(
        os.path.join(out_dir, "{}_norm_log.txt".format(prefix)),
        float_format="%0.6f",
        sep="\t",
        header=True,
        columns=sorted_columns,
        index=True,
    )


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
        "-m",
        "--mirbase_version",
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
        "-c",
        "--config_file",
        help="Custom configuration file",
    )

    parser.add_argument(
        "-o",
        "--out_dir",
        help="output directory path for expression matrices",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
