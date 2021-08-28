import argparse
import os
import logging
import pandas

from pathlib import Path

from mirna_profiling.library_stat.graph_stat import generate_project_graphs
from mirna_profiling.utility.csv_util import combine_csv

EXP_FILES = [
    "{}{}".format(x, y)
    for x in ["expn_matrix", "expn_matrix_mimat"]
    for y in [".txt", "_norm.txt", "_norm_log.txt"]
]


def main():

    args = parse_args()
    if args.out_dir:
        out_dir = args.out_dir
    else:
        out_dir = args.proj_dir

    project_summary(args.proj_dir, out_dir, proj_name=args.name)


def project_summary(proj_dir, out_dir, proj_name=None):
    os.makedirs(out_dir, 0o775, exist_ok=True)
    logging.info("Collating alignment stats")
    combine_csv(
        collate_sample_files(proj_dir, "alignment_stats.csv"),
        axis=0,
        index_col="Library",
        header=0,
        out_csv=os.path.join(out_dir, "alignment_stats.csv"),
    )
    logging.info("Collating expression matrix files")
    for exp_file in EXP_FILES:
        combine_csv(
            collate_sample_files(proj_dir, exp_file),
            axis=1,
            index_col="Gene",
            sep="\t",
            header=0,
            out_csv=os.path.join(out_dir, exp_file),
            sort_column=True,
            sort_index=True,
        )

    logging.info("Generating graphs")
    graph_dir = os.path.join(out_dir, "graphs")
    os.makedirs(graph_dir, 0o775, exist_ok=True)
    generate_project_graphs(proj_dir, out_dir=graph_dir, proj_name=proj_name)


def collate_sample_files(proj_dir, fname):

    csv_files = []
    for feature_dir in Path(proj_dir).rglob("*_features"):
        sample = str(os.path.basename(feature_dir)).replace("_features", "")
        sample_dir = os.path.dirname(feature_dir)
        alignment_stat_csv = os.path.join(sample_dir, fname)
        if os.path.exists(alignment_stat_csv):
            csv_files.append(alignment_stat_csv)
        else:
            logging.warning(
                "Missing {} for sample {} in {}".format(fname, sample, sample_dir)
            )

    return csv_files


def parse_args():
    parser = argparse.ArgumentParser(description="Generate graphs for miRNA project")
    parser.add_argument(
        "-p", "--proj_dir", help="miRNA project directory", required=True
    )
    parser.add_argument("-n", "--name", help="name to use in graphs", required=True)
    parser.add_argument("-o", "--out_dir", help="Save graphs in output directory")
    return parser.parse_args()


if __name__ == "__main__":
    main()
