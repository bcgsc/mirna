#!/usr/bin/env python

import argparse
import glob
import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy
import pandas
from scipy import interpolate
from statsmodels.nonparametric.smoothers_lowess import lowess

from mirna_profiling.utility.logging_util import init_logger

NORM_FACTOR = 1000000

logger = init_logger()


def main():
    args = parse_args()
    generate_project_graphs(
        os.path.abspath(args.proj_dir),
        plot_type=args.plot_type,
        out_dir=args.out_dir,
        proj_name=args.name,
    )


def generate_project_graphs(proj_dir, plot_type=None, out_dir=None, proj_name=None):
    """
    Generate summary alignment_stats, expression matrix and plots for all samples in proj_dir

    Args:
        proj_dir (str): output directory
        plot_type (str): saturation, tag or adapter
        out_dir (str): output directory

    Returns:

    """
    if not out_dir:
        out_dir = proj_dir
    tag_dir = os.path.join(out_dir, "tag")
    adapter_dir = os.path.join(out_dir, "adapter")
    for mdir in (out_dir, tag_dir, adapter_dir):
        os.makedirs(mdir, 0o775, exist_ok=True)

    if not proj_name:
        proj_name = str(os.path.basename(os.path.realpath(proj_dir)))
    alignment_stat = glob.glob("{}/alignment_stats.csv".format(proj_dir))
    if alignment_stat and (not plot_type or plot_type == "saturation"):
        plot_saturation(alignment_stat[0], proj_name, out_dir=out_dir)

    for feature_dir in Path(proj_dir).rglob("*_features"):
        sample_name = str(os.path.basename(feature_dir)).replace("_features", "")
        sample_dir = os.path.dirname(feature_dir)
        if not plot_type or plot_type == "adapter":
            adapter_file = glob.glob(
                "{}/{}*adapter*report*".format(sample_dir, sample_name)
            )
            if adapter_file:
                plot_adapter(adapter_file[0], libname=sample_name, out_dir=adapter_dir)
        if not plot_type or plot_type == "tag":
            taglength_files = Path(feature_dir).glob("*_taglengths.csv")
            for taglen in taglength_files:
                plot_taglengths(taglen, libname=sample_name, out_dir=tag_dir)


def plot_adapter(adapter_file, libname=None, out_dir=None):
    if not libname:
        libname = os.path.basename(adapter_file).split("_")[0]
    if not out_dir:
        out_dir = os.getcwd()

    df = pandas.read_csv(adapter_file, sep=" ", header=None, index_col=0).sort_index()

    ax = df.plot.bar(
        stacked=False,
        xlabel="tag length (bp)",
        ylabel="number of reads",
        title="{} - Position Along Index Trimmed Read With 3' Adapter".format(libname),
        rot=90,
        width=0.85,
        legend=False,
    )
    ax.title.set_size(10)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "{}_adapter_report.jpeg".format(libname)))
    plt.close()


def plot_saturation(alignment_stat, proj_name, out_dir=None):
    if not out_dir:
        out_dir = os.getcwd()

    data = pandas.read_csv(alignment_stat)
    data["Total miRNA"] = data["Total miRNA"] / NORM_FACTOR

    plt.scatter(
        data=data,
        x="Total miRNA",
        y="miRNA Species",
        label=">= 1x coverage",
        facecolors="none",
        edgecolors="blue",
    )
    plt.scatter(
        data=data,
        x="Total miRNA",
        y="miRNA Species Covered by >= 10 Reads",
        label=">= 10x coverage",
        facecolors="none",
        edgecolors="red",
    )
    sat1_xnew, sat1_ynew = smoothing(data["Total miRNA"], data["miRNA Species"])
    sat10_xnew, sat10_ynew = smoothing(
        data["Total miRNA"], data["miRNA Species Covered by >= 10 Reads"]
    )
    if sat1_xnew is not None:
        plt.plot(
            sat1_xnew,
            sat1_ynew,
            color="blue",
        )

    if sat10_xnew is not None:
        plt.plot(
            sat10_xnew,
            sat10_ynew,
            color="red",
        )

    plt.legend(loc="lower right")
    plt.xlabel("# reads aligned to miRNAs (Millions)")
    plt.ylabel("# miRNA species")
    plt.title("miRNA Saturation in {}".format(proj_name))
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    plt.savefig(os.path.join(out_dir, "{}_saturation.jpeg".format(proj_name)))
    plt.close()


def smoothing(x, y, lowess_frac=0.6, sparse=100, min_samples=4):
    xsmoothed = None
    ysmoothed = None
    try:
        # two LOWESS curves
        if len(x) >= min_samples:
            logger.info("LOWESS smoothing of dense data")
            lowessed = lowess(y, x, frac=lowess_frac)
            xsmoothed = lowessed[:, 0]
            ysmoothed = lowessed[:, 1]
        # if len(x) >= min_samples:
        #     logger.info("BSpline smoothing of sparse data")
        #     x_array = x.sort_values()
        #     # 300 represents number of points to make between .min and .max
        #     xsmoothed = numpy.linspace(x_array.min(), x_array.max(), 300)
        #     # generate BSpline object
        #     spl = interpolate.make_interp_spline(x_array, y, k=3)
        #     ysmoothed = spl(xsmoothed)
        else:
            logging.warning(
                "No smoothing available for fewer than {} samples".format(min_samples)
            )
    except Exception as e:
        logging.warning("Smoothing failure: {}".format(e))
        xsmoothed = None
        ysmoothed = None
    return xsmoothed, ysmoothed


def plot_taglengths(taglen_csv, libname, out_dir=None):
    prefix = os.path.basename(taglen_csv).split("_")[0]
    if not out_dir:
        out_dir = os.getcwd()
    data = pandas.read_csv(taglen_csv)
    if data.empty:
        logger.warning("{} is empty".format(taglen_csv))
        return
    data.fillna(0, inplace=True)
    data["taglen"].astype(int)
    data.set_index("taglen", inplace=True)
    rna = (
        data["snoRNA"]
        + data["tRNA"]
        + data["rRNA"]
        + data["snRNA"]
        + data["scRNA"]
        + data["srpRNA"]
        + data["rmsk_RNA"]
        + data["No_CDS"]
        if "No_CDS" in data.columns
        else data["Non-Coding_Exon"]
    )
    coding_gene = (
        (data["p5_UTR"] if "p5_UTR" in data.columns else data["5_UTR"])
        + (data["p3_UTR"] if "p3_UTR" in data else data["3_UTR"])
        + data["Coding_Exon"]
        + data["Intron"]
    )
    repeat = (
        data["LINE"]
        + data["SINE"]
        + data["LTR"]
        + data["Satellite"]
        + data["rmsk_DNA"]
        + data["rmsk_Low_complexity"]
        + data["rmsk_Simple_repeat"]
        + data["rmsk_Other"]
        + data["rmsk_Unknown"]
    )

    df = pandas.concat(
        [
            data["mature"],
            data["star"],
            data["unannotated"],
            data["crossmapped"],
            data["stemloop"],
            data["precursor"],
            rna,
            coding_gene,
            repeat,
            data["Unknown"],
        ],
        axis=1,
    )

    mlabels = [
        "miRNA, mature strand",
        "miRNA, star strand",
        "miRNA, unannotated in miRBase",
        "miRNA, crossmapped",
        "miRNA, stemloop",
        "miRNA, precursor",
        "Other RNAs",
        "Coding Genes",
        "Repeats",
        "Unannotated",
    ]

    mcolors = [
        "blue",
        "royalblue",
        "darkblue",
        "slateblue",
        "purple",
        "darkslateblue",
        "red",
        "burlywood",
        "chartreuse",
        "grey",
    ]
    ax = df.plot.bar(
        stacked=True,
        xlabel="tag length (bp)",
        ylabel="% of all reads aligning to <3 positions",
        title="{} - Percentage of Aligned Tags At Each Tag Length With Annotation".format(
            prefix
        ),
        rot=0,
        color=mcolors,
        width=0.85,
    )
    ax.legend(loc="best", fontsize="x-small", labels=mlabels)
    plt.savefig(os.path.join(out_dir, "{}_{}_taglengths.jpeg".format(libname, prefix)))
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser(description="Generate graphs for miRNA project")
    parser.add_argument(
        "-p", "--proj_dir", help="miRNA project directory", required=True
    )
    parser.add_argument("-o", "--out_dir", help="Save graphs in output directory")
    parser.add_argument(
        "-t",
        "--plot_type",
        help="Type of plots to generate. Default=%(default)s",
        choices=["saturation", "tag", "adapter"],
    )
    parser.add_argument("-n", "--name", help="name of project to show in plots")
    return parser.parse_args()


if __name__ == "__main__":
    main()
