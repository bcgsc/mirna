#!/usr/bin/env python

import argparse
import glob
import gzip
import logging
import os
import re
from pathlib import Path

import pandas
import pysam
import yaml

from mirna_profiling.annotation.annotate_bam import (
    parse_cigar_string,
    is_annotated,
    annotate_sam,
)
from mirna_profiling.annotation.feature_db import FeatureDatabase
from mirna_profiling.annotation.mirbase import Mirbase
from mirna_profiling.annotation.ucsc import UCSC
from mirna_profiling.library_stat.adapter import get_adapter_stat
from mirna_profiling.library_stat.tcga import write_tcga_files
from mirna_profiling.utility.logging_util import init_logger

annotation_tags = ["XC", "XI", "XD"]
default_tag_type = "Z"


flag_columns = [
    "aligned",
    "filtered",
    "chastity_failed",
    "overfiltered",
    "softclipped",
    "crossmapped",
    "mirna_mapped",
]

stat_columns = flag_columns + [
    "feature_type",
    "feature_name",
    "mirna_species",
    "taglen",
]

logger = init_logger("{}.log".format(__name__))


class AlignmentStat(object):
    """
    A class to calculate the alignment stats of a sample/bam or samples/bams in a project (directory)

    TODO: Reduce memory usage
    TODO: Standardize library name extraction from bam name, first word or token split with . and _?

    """

    def __init__(self, config_file=None):
        """
        Load config file and initialize object

        # TODO: Validate annotation config

        """

        if not config_file:
            config_file = os.path.join(
                os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
                "configuration",
                "annotation_config.yaml",
            )

        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)

        self.feature_priority = self.config["Feature Annotation"]["feature_priority"]
        self.feature_type_csv_header = self.config["Feature Annotation"].get(
            "feature_type_csv_header"
        )

        self.MAX_NM = self.config["Alignment Stat"].get("MAX_NM")
        self.MAX_X0 = self.config["Alignment Stat"].get("MAX_X0")
        self.MAX_XA = self.config["Alignment Stat"].get("MAX_XA")

    def get_read_stats(self, inbam, mirbase, feature_dir):
        """
        Generate raw read stats, including isoforms, on bam and bed files for peak discovery

        Args:
            inbam (str): input annotated bam
            bed_dir (str): directory for bed files
            mirbase (Mirbase): mirbase for isoform nomenclature

        Returns:
            (tuple): read stats and isoforms

        """

        logger.info("Started generating raw read stats on b/sam file: {}".format(inbam))

        read_stats = []
        isoforms = {}
        bedfile_objs = {}
        samfile = pysam.AlignmentFile(inbam, "r")
        progress = 0
        bed_dir = os.path.join(feature_dir, "bed")
        os.makedirs(bed_dir, 0o775, exist_ok=True)
        for read in samfile.fetch(until_eof=True):
            read_stat, alns, mirna_indices = self.get_read_stat(read)
            read_stats.append(read_stat)
            if read_stat.get("filtered") and not read_stat.get("softclipped"):
                if read_stat.get("mirna_mapped"):
                    mirnas = [alns[x] for x in mirna_indices]
                    self.record_isoform(isoforms, mirnas, read.query_alignment_sequence)
                    if read_stat.get("crossmapped"):
                        self.write_bedfile(bedfile_objs, bed_dir, mirnas)
                if not read_stat.get("crossmapped"):
                    self.write_bedfile(bedfile_objs, bed_dir, (alns[0],))
            progress += 1
            if progress % 100000 == 0:
                logging.info(f"Processed {progress:,} reads")

        for key, fobj in bedfile_objs.items():
            fobj.close()
        samfile.close()

        self.update_isoform_nomenclature(isoforms, mirbase)
        self.write_isoform(isoforms, os.path.join(feature_dir, "isoforms.txt"))

        logger.info(
            "Finished generating raw read stats on b/sam file: {}".format(inbam)
        )

        return read_stats

    def get_sample_stat(
        self,
        libname,
        inbam,
        genome,
        mirbase,
        adapter_file=None,
        out_dir=None,
    ):
        """
        Generate alignment stat and feature expression files in sample output dir

        Args:
            libname (str): library name
            inbam (str): input annotated bam
            genome (str): UCSC genome name, eg. hg38
            mirbase (Mirbase): an instance of Mirbase class
            out_dir (str): output directory
            adapter_file (str): adapter report file

        Returns:

        """
        if not is_annotated(inbam):
            logger.warning(
                "{} has not been annotated. Please annotate the bam before running alignment_stat".format(
                    inbam
                )
            )
            return None

        logger.info("Started generating alignment stat on b/sam file: {}".format(inbam))
        if not out_dir:
            out_dir = Path(inbam).parent

        if not adapter_file:
            adapter_file = list(
                Path(inbam).parent.glob("{}*adapter*report*".format(libname))
            )
            if adapter_file:
                adapter_file = adapter_file[0]

        feature_dir = os.path.join(out_dir, "{}_features".format(libname))

        read_stats = self.get_read_stats(inbam, mirbase, feature_dir)

        df = pandas.DataFrame(read_stats)
        crossmapped_features = (
            df.loc[df["crossmapped"] == True, "feature_name"]
            .apply(lambda x: x.split(";"))
            .explode("feature_name")
            .unique()
        )
        df.loc[df["feature_name"].isin(crossmapped_features), "crossmapped"] = True
        df["chastity_failed"] = df.apply(
            lambda x: (x["chastity_passed"] == False) & (x["filtered"] == True), axis=1
        )
        total_reads = len(df)
        min_taglen = df["taglen"].dropna().astype(int).min()
        bstat = df[flag_columns].agg("sum").astype(int).to_dict()
        bam_stat = {
            "Library": libname,
            "Total Reads >= {}bp".format(min_taglen): total_reads,
        }
        if adapter_file:
            logger.info("Adapter file is present: {}".format(adapter_file))
            adapter_stat = get_adapter_stat(adapter_file)
            bam_stat.update(adapter_stat)
        else:
            logger.info("Adapter file is not present with {}".format(inbam))
        bam_stat.update(
            {
                "Aligned Reads Post-Filter": bstat["filtered"],
                "% Aligned Reads": "{:0.2f}%".format(
                    bstat["filtered"] * 100 / total_reads
                ),
                "Unaligned Reads": total_reads - bstat["aligned"],
                "% Unaligned Reads": "{:0.2f}%".format(
                    (total_reads - bstat["aligned"]) * 100 / total_reads
                ),
                "Filtered Reads without XA": bstat["overfiltered"],
                "Softclipped Reads": bstat["softclipped"],
                "Chastity Failed Reads Post-Filter": bstat["chastity_failed"],
                "miRNA Species": 0,
                "miRNA Species Covered by >= 10 Reads": 0,
                "Total miRNA": bstat["mirna_mapped"],
                "Crossmapped miRNA": bstat["crossmapped"],
            }
        )

        df_mirna_species = (
            df.loc[df["mirna_mapped"] == True, ["mirna_species", "mirna_mapped"]]
            .explode("mirna_species")
            .groupby("mirna_species")
            .count()
        )

        bam_stat["miRNA Species"] = df_mirna_species.shape[0]
        bam_stat["miRNA Species Covered by >= 10 Reads"] = df_mirna_species[
            df_mirna_species["mirna_mapped"] >= 10
        ].shape[0]

        for feature_type in self.feature_priority + ["Unknown"]:
            display_feature_type = (
                self.feature_type_csv_header.get(feature_type)
                if self.feature_type_csv_header
                and feature_type in self.feature_type_csv_header
                else feature_type
            )
            bam_stat[display_feature_type] = df.loc[
                (df["feature_type"] == feature_type) & (df["crossmapped"] != True),
                "filtered",
            ].count()

        bam_stat["% Total miRNA"] = "{:0.2f}%".format(
            bam_stat["Total miRNA"] * 100 / bam_stat["Aligned Reads Post-Filter"]
        )

        bam_stat["% Crossmapped miRNA"] = "{:0.2f}%".format(
            bam_stat["Crossmapped miRNA"] * 100 / bam_stat["Aligned Reads Post-Filter"]
        )

        for feature_type in self.feature_priority + ["Unknown"]:
            display_feature_type = (
                self.feature_type_csv_header.get(feature_type)
                if self.feature_type_csv_header
                and feature_type in self.feature_type_csv_header
                else feature_type
            )
            bam_stat["% {}".format(display_feature_type)] = "{:0.2f}%".format(
                bam_stat[display_feature_type]
                * 100
                / bam_stat["Aligned Reads Post-Filter"]
            )

        self.write_mirna_files(
            df_mirna_species,
            os.path.join(feature_dir, "mirna_species.txt"),
            group_by="mirna_species",
            norm_factor=bam_stat["Total miRNA"],
        )

        self.write_mirna_files(
            df.loc[df["crossmapped"] == True, ["mirna_mapped", "feature_name"]],
            os.path.join(feature_dir, "crossmapped.txt"),
            norm_factor=bam_stat["Total miRNA"],
        )

        self.write_mirna_files(
            df.loc[
                (df["mirna_mapped"] == True) & (df["crossmapped"] != True),
                ["mirna_mapped", "feature_name"],
            ],
            os.path.join(feature_dir, "miRNA.txt"),
            norm_factor=bam_stat["Total miRNA"],
        )

        self.write_taglengths(
            df,
            "filtered",
            os.path.join(feature_dir, "filtered_taglengths.csv"),
            norm_factor=bam_stat["Aligned Reads Post-Filter"],
        )

        self.write_taglengths(
            df,
            "chastity_failed",
            os.path.join(feature_dir, "chastity_taglengths.csv"),
            norm_factor=bam_stat["Chastity Failed Reads Post-Filter"],
        )
        self.write_taglengths(
            df,
            "softclipped",
            os.path.join(feature_dir, "softclip_taglengths.csv"),
            norm_factor=bam_stat["Softclipped Reads"],
        )

        self.write_feature_species(df, feature_dir)

        pandas.DataFrame([bam_stat]).to_csv(
            os.path.join(out_dir, "alignment_stats.csv"),
            float_format="%0.2f",
            header=True,
            index=False,
        )

        write_tcga_files(out_dir, genome, mirbase)

        logger.info(
            "Finished generating alignment stat on b/sam file: {}".format(inbam)
        )
        return bam_stat

    @staticmethod
    def record_isoform(isoforms, mirnas, read_sequence):
        mirna_count = len(mirnas)
        for aln in mirnas:
            coord = "{} {} {} {}".format(
                aln.get("chr"),
                aln.get("pos"),
                aln.get("pos") + aln.get("taglen") - 1,
                aln.get("strand"),
            )
            mirna_species = get_name(aln)
            feature_key = get_feature_key(aln)
            if mirna_species not in isoforms:
                isoforms[mirna_species] = {}
            if coord not in isoforms.get(mirna_species):
                isoforms[mirna_species][coord] = {}
                isoforms[mirna_species][coord]["count"] = 1
                isoforms[mirna_species][coord]["crossmapped"] = (
                    1 if mirna_count > 1 else 0
                )
                isoforms[mirna_species][coord]["seq"] = read_sequence
            else:
                isoforms[mirna_species][coord]["count"] += 1
                if mirna_count > 1:
                    isoforms[mirna_species][coord]["crossmapped"] = 1
            isoforms[mirna_species][coord]["annot"] = ",".join(
                feature_key.split(",")[1:]
            )

    @staticmethod
    def update_isoform_nomenclature(isoforms, mirbase):
        """

        Args:
            isoforms (dict):
            mirbase (Mirbase):

        Returns:
            isoforms(dict): isoforms with updated nomenclature
        """

        for mirna_species, mdict in isoforms.items():
            for coord in mdict.keys():
                (chr, isomir_start, isomir_end, strand) = coord.split()
                isomir_start = int(isomir_start)
                isomir_end = int(isomir_end)
                if "mature" in mdict[coord]["annot"]:
                    mature_name = mdict[coord]["annot"].split(",")[1]
                    mature_mirnas = [
                        x
                        for x in mirbase.get_mature_mirnas_by_id(mature_name)
                        if x.get("start") <= isomir_end and x.get("end") >= isomir_start
                    ]
                    if mature_mirnas:
                        mature_mirna = mature_mirnas[0]
                    else:
                        mature_mirna = None
                    if not mature_mirna:
                        logger.warning(
                            "mature miRNA {} not in mirbase".format(mature_name)
                        )
                        mdict[coord]["isomir_name"] = "NA"
                    else:
                        mdict[coord]["isomir_name"] = "|".join(
                            [
                                mature_mirna.get("mature_name"),
                                prefix_sign(isomir_start - mature_mirna.get("start")),
                                prefix_sign(isomir_end - mature_mirna.get("end")),
                                "",
                            ]
                        )
                else:
                    mdict[coord]["isomir_name"] = "NA"

        return isoforms

    @staticmethod
    def write_isoform(isoforms, fpath):
        with open(fpath, "w") as f:
            for mirna_species, mdict in isoforms.items():
                for coord in mdict.keys():
                    f.write(
                        "{} {} {} {} {} {} {}\n".format(
                            mirna_species,
                            coord,
                            mdict[coord]["seq"],
                            mdict[coord]["count"],
                            mdict[coord]["crossmapped"],
                            mdict[coord]["annot"],
                            mdict[coord]["isomir_name"],
                        )
                    )

    @staticmethod
    def write_bedfile(file_objs, out_dir, alns):
        """

        Args:
            file_objs (dict)): bed file objects
            alns (iterable): alignments

        Returns:

        """
        for aln in alns:
            fname = "{}_{}".format(aln.get("chr"), aln.get("strand"))
            if fname not in file_objs:
                outfile = os.path.join(out_dir, "{}.txt.gz".format(fname))
                file_objs[fname] = gzip.open(outfile, "wt")
            file_objs[fname].write(
                "{}\t{}\t{}\n".format(
                    aln.get("chr"),
                    aln.get("pos") - 1,
                    aln.get("pos") + aln.get("taglen") - 1,
                )
            )

    def get_read_stat(self, read):
        """
        Get stat for a read (psam.AlignmentSegment)

        Args:
            read (pysam.AlignedSegment):

        Returns:
            read_stat, alns, mirna_indices (tuple):
                read_stat(dict):
                    "aligned"
                    "filtered"
                    "chastity_passed"
                    "overfiltered"
                    "softclipped"
                    "crossmapped"
                    "mirna_mapped"
                    "mirna_species"
                    "feature_type"
                    "feature_name"
                    "taglen"
        """
        read_stat = dict.fromkeys(flag_columns, False)
        read_stat.update(dict.fromkeys(stat_columns))

        if not read.reference_name or read.reference_name == "*":
            read_stat["aligned"] = False
            return read_stat, None, None
        else:
            read_stat["aligned"] = True

        nm = get_tag_value(read, "NM", "i")
        x0 = get_tag_value(read, "X0", "i")
        if nm > self.MAX_NM or x0 > self.MAX_X0:
            read_stat["filtered"] = False
            return read_stat, None, None
        else:
            read_stat["filtered"] = True

        x1 = get_tag_value(read, "X1", "i")

        read_stat["overfiltered"] = (
            True if x0 <= self.MAX_X0 and x0 + x1 > self.MAX_XA else False
        )

        read_stat["chastity_passed"] = False if (read.flag & 0x0200) == 0x0200 else True

        alns = self.get_annotated_alignments(read)
        annot_key = get_annot_key(alns[0])

        read_stat["taglen"] = alns[0].get("taglen")
        read_stat["feature_type"] = (
            annot_key if annot_key in self.feature_priority else "Unknown"
        )

        if alns[0].get("is_softclipped") == True:
            read_stat["softclipped"] = True
            return read_stat, None, None
        else:
            read_stat["softclipped"] = False

        feature_key = get_feature_key(alns[0])
        read_stat["feature_name"] = feature_key
        mirnas = self.get_mirna_alignments(alns)

        if mirnas:
            read_stat["mirna_mapped"] = True
            if len(mirnas) > 1:
                read_stat["crossmapped"] = True
                read_stat["feature_type"] = "crossmapped"
                read_stat["mirna_species"] = [get_name(alns[x]) for x in mirnas]
                long_key = ";".join(sorted(get_feature_key(alns[x]) for x in mirnas))
                read_stat["feature_name"] = long_key
            else:
                read_stat["crossmapped"] = False
                read_stat["mirna_species"] = [get_name(alns[mirnas[0]])]
                read_stat["feature_name"] = get_feature_key(alns[mirnas[0]])

        else:
            read_stat["mirna_mapped"] = False

        return read_stat, alns, mirnas

    def get_annotated_alignments(self, read):
        """
        Get alignments of a read, sorted by annotated feature priority

        Args:
            read (pysam.AlignedSegment)

        Returns:
            alns (list):
        """

        x0 = get_tag_value(read, "X0", "i")
        reference_start = read.reference_start + 1
        aln = "{},{},{},{}".format(
            read.reference_name,
            -1 * reference_start if (read.flag & 0x0010) == 0x0010 else reference_start,
            read.cigarstring,
            get_tag_value(read, "NM", "i"),
        )

        if read.has_tag("XA"):
            xa_tag = "{};{};".format(aln, get_tag_value(read, "XA"))
        else:
            xa_tag = "{};".format(aln)

        read.set_tag("XA", xa_tag)

        use_tags = ["XA"] + annotation_tags

        xs = (
            dict(zip(use_tags, x))
            for x in zip(*(get_tag_value(read, x).split(";") for x in use_tags))
        )

        alns = []
        num_aln = 0
        for x in xs:
            if not x["XA"]:
                continue
            xa_items = x["XA"].split(",")
            if int(xa_items[3]) > self.MAX_NM:
                continue
            x["chr"] = xa_items[0]
            x["pos"] = int(xa_items[1].replace("+", ""))
            if x["pos"] < 0:
                x["strand"] = "-"
                x["pos"] *= -1
            else:
                x["strand"] = "+"
            x["is_softclipped"] = True if "S" in xa_items[2] else False
            x["taglen"] = parse_cigar_string(xa_items[2])
            alns.append(x)
            num_aln += 1
            if num_aln == x0:
                break

        sorted_alns = sorted(
            alns,
            key=lambda x: self.feature_priority.index(get_annot_key(x))
            if get_annot_key(x) in self.feature_priority
            else 1000000,
        )

        return sorted_alns

    @staticmethod
    def get_mirna_alignments(alns):
        """
        Get aligned miRNAs from a list of annotated alignments

        Args:
            alns (list):

        Returns:
            mirnas (list): indices of miRNA alignments in input alignments
        """

        mirnas = []
        mirna_clans = {}
        for idx in range(0, len(alns)):
            mirna_clan = None
            m = re.match(r"[^-]+-[^-]+-[^-]+", get_name(alns[idx]))
            if m is not None:
                mirna_clan = m.group(0)
            if (
                alns[idx].get("XC") == "miRNA"
                and alns[idx].get("is_softclipped") == False
                and mirna_clan not in mirna_clans
            ):
                mirna_clans[mirna_clan] = 1
                mirnas.append(idx)

        return mirnas

    @staticmethod
    def write_mirna_files(df, fpath, group_by="feature_name", norm_factor=None):
        if df.empty:
            logging.warning("Dataset is empty. Cannot generate file: {}".format(fpath))
            return

        df_mirna = (
            df.groupby(group_by)
            .sum()
            .astype(int)
            .fillna(0)
            .sort_values("mirna_mapped", ascending=False)
        )
        if norm_factor:
            df_mirna["% of Total miRNA"] = df_mirna["mirna_mapped"] * 100 / norm_factor
        df_mirna.to_csv(
            fpath, float_format="%0.2f%%", sep=" ", index=True, header=False
        )

    def write_taglengths(self, df, column_name, fpath, norm_factor=None):

        df_filtered = (
            df.loc[df[column_name] == True, [column_name, "feature_type", "taglen"]]
            .astype({"taglen": "int64"})
            .groupby(["feature_type", "taglen"], as_index=False)
            .sum()
        )

        if df_filtered.shape[0] > 0:
            df_filtered = df_filtered.pivot("taglen", "feature_type").fillna(0)
            df_filtered.columns = (x[1] for x in df_filtered.columns)
            for feature_type in ["crossmapped"] + self.feature_priority:
                if feature_type not in df_filtered.columns:
                    df_filtered[feature_type] = 0
        else:
            df_filtered = pandas.DataFrame(
                {x: [] for x in ["taglen", "crossmapped"] + self.feature_priority}
            )
            df_filtered.set_index("taglen", inplace=True)

        df_filtered.columns = (x.replace(" ", "_") for x in df_filtered.columns)
        if norm_factor:
            df_filtered = df_filtered.applymap(lambda x: x * 100 / norm_factor)
        df_filtered.to_csv(fpath, header=True)

    def write_feature_species(self, df, out_dir):
        for feature_type in self.feature_priority:
            df_features = (
                df.loc[
                    (df["feature_type"] == feature_type) & (df["mirna_mapped"] != True)
                ]
                .groupby(["feature_name"])["filtered"]
                .sum()
                .astype(int)
                .fillna(0)
            ).sort_values(ascending=False)
            if len(df_features) > 0:
                df_features.to_csv(
                    os.path.join(
                        out_dir, "{}.txt".format(feature_type.replace(" ", "_"))
                    ),
                    sep=" ",
                    header=False,
                )


def prefix_sign(number):

    if number > 0:
        signed_number = "+{}".format(number)
    else:
        signed_number = "{}".format(number)

    return signed_number


def get_annot_key(x):
    return x.get("XD").split(",")[0] if x.get("XD") else x.get("XC")


def get_feature_key(x):
    if x.get("XD"):
        key = "{},{}".format(x.get("XI"), x.get("XD"))
    else:
        key = x.get("XI")

    return key


def get_name(x):
    return x.get("XI")


def get_tag_value(read, tag_name, tag_code="Z"):
    if tag_code == "i":
        tag_value = int(read.get_tag(tag_name)) if read.has_tag(tag_name) else 0
    else:
        tag_value = read.get_tag(tag_name) if read.has_tag(tag_name) else ""

    return tag_value


def get_project_stat(
    proj_dir,
    ucsc_version,
    mirbase_version,
    species_code,
    out_dir=None,
    config_file=None,
    no_subdir=False,
    keep_bam=False,
):
    """
    Get alignment stats for all bams in a project directory

    Args:
        proj_dir (str) : path of bam file or project directory containing bam files
        ucsc_version (str) : UCSC genome or databae version, e.g. hg38
        mirbase_version (str): miRBase version, e.g. mirna_21
        species_code (str): miRBase species coce, e.g. hsa
        out_dir (str): output directory path
        config_file (str): custom config file
        keep_bam (bool): keep annotated b/sam  file

    Returns:

    """
    if os.path.isdir(proj_dir):
        bams = [str(x) for x in Path(proj_dir).rglob("*.[bs]am")]
    else:
        bams = [proj_dir]

    aln_stat = AlignmentStat(config_file=config_file)
    proj_stats = []
    feature_db = UCSC(ucsc_version, config_file=config_file)
    mirbase = Mirbase(mirbase_version, species_code, config_file=config_file)
    feature_db.add_features(mirbase.get_features())

    for inbam in bams:
        fname = str(os.path.basename(inbam).split(".")[0])
        tokens = fname.split("_")
        libname = tokens[0]

        if len(tokens) > 1:
            m = re.match("^([ATCG]+)$", tokens[1])
            if m:
                libname += "_" + m.group(1)

        if not out_dir:
            sample_dir = Path(inbam).parent
        elif no_subdir and len(bams) == 1:
            sample_dir = out_dir
        else:
            sample_dir = os.path.join(out_dir, libname)

        logger.info("Processing b/sam file: {}".format(inbam))
        if not is_annotated(inbam):
            logger.info("Missing annotation in b/sam file: {}".format(inbam))
            do_annot = True
            os.makedirs(sample_dir, 0o775, exist_ok=True)
            annotated_bam = os.path.join(sample_dir, os.path.basename(inbam) + ".annot")
            annotate_sam(feature_db, inbam, annotated_bam, config_file)
        else:
            do_annot = False
            logger.info("B/sam has been annotated already: {}".format(inbam))
            annotated_bam = inbam

        adapter_file = list(
            Path(inbam).parent.glob("{}*adapter*report*".format(libname))
        )
        if adapter_file:
            adapter_file = adapter_file[0]
        else:
            adapter_file = None

        sample_stat = aln_stat.get_sample_stat(
            libname,
            annotated_bam,
            ucsc_version,
            mirbase,
            adapter_file=adapter_file,
            out_dir=sample_dir,
        )
        if sample_stat:
            proj_stats.append(sample_stat)

        if not keep_bam and do_annot and os.path.exists(annotated_bam):
            os.remove(annotated_bam)

        if adapter_file:
            symlink_file(adapter_file, sample_dir, "{}_adapter.report".format(libname))

    pandas.DataFrame(proj_stats).to_csv(
        os.path.join(out_dir, "alignment_stats.csv"),
        float_format="%0.2f",
        header=True,
        index=False,
    )

    status_file = os.path.join(
        out_dir,
        "alignment_stat.{}".format(
            "success" if len(proj_stats) == len(bams) else "failure"
        ),
    )
    Path(status_file).touch()

    return proj_stats


def symlink_file(infile, target_dir, link=None):
    """
    Create symlink for infile in target_dir

    Args:
        infile: source file path
        link: symlink name if different from source file name
        target_dir: target directory to create the symlink

    """
    do_link = False
    infile = Path(infile).resolve()
    sname = os.path.basename(infile)
    target_dir = Path(target_dir).resolve()
    if infile.parent.samefile(target_dir):
        if link and link != sname:
            do_link = True
    else:
        if not link:
            link = sname
        do_link = True
    if do_link:
        target_file = Path(os.path.join(target_dir, link))
        if not target_file.exists():
            target_file.symlink_to(infile)


def main():
    args = parse_args()
    # print(timeit.timeit(get_project_stat, number=1))
    if args.out_dir:
        out_dir = args.out_dir
    elif os.path.isdir(args.bam):
        out_dir = args.bam
    else:
        out_dir = os.getcwd()

    get_project_stat(
        args.bam,
        args.ucsc,
        args.mirbase,
        args.species,
        out_dir=out_dir,
        config_file=args.config_file,
        no_subdir=args.nd,
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Generate alignment stat")
    parser.add_argument(
        "-u",
        "--ucsc",
        help="UCSC genome/database name, e,g. hg38, hg19, mm10",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--mirbase",
        help="miRBase MySQL database name, e.g. mirna_21, mirna_16",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--species",
        help="miRBase species code, e.g., hsa for human, mmu for mouse",
        required=True,
    )
    group_annot = parser.add_mutually_exclusive_group(required=True)
    group_annot.add_argument(
        "-p",
        "--proj_dir",
        dest="bam",
        help="Project directory path with annotated b/sam files",
    )
    group_annot.add_argument(
        "-b",
        "--bam",
        dest="bam",
        help="annotated b/sam file",
    )
    parser.add_argument(
        "-c",
        "--config_file",
        help="Custom config file",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        help="Output directory",
    )
    parser.add_argument(
        "--nd",
        help="Do not create sample directories under out_dir",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-k",
        "--keep_bam",
        help="Keep annotated bam file",
        default=False,
        action="store_true",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
