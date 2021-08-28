#!/usr/bin/env python

import json
import logging
import os
import timeit
from string import Template

# Password authentication problem with mysql-connector-python for MySQL 5.5 and older
# Version mismatch ("Bad handshake") with mysql-connector-python for MySQL 5.0
# import mysql.connector as sql
import pandas
import pymysql as sql

from mirna_profiling.annotation.feature_db import FeatureDatabase, load_csv
from mirna_profiling.utility.gff_util import gff_to_records, GFF_COLUMNS

resource_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "resource")


class Mirbase(FeatureDatabase):
    """
    A class that acts as an interface to miRBase, and retrieves features of interest from the miRBase MySQL databases

    """

    mirna_attrs = [
        "mirna_id",
        "chr",
        "start",
        "end",
        "strand",
    ]

    mature_mirna_attrs = ["mature_from", "mature_to", "mature_name", "mature_acc"]

    def __init__(
        self,
        database,
        species_code,
        user=None,
        host=None,
        password=None,
        **kwarg,
    ):
        """
        Initialize config parameters and data resources

        Args:
            database (str): miRBase release version or database name following the pattern "mirna_<release version>"
            species_code (str): miRBase species code, e.g. hsa for human
            **kwarg:
        """
        super().__init__(**kwarg)

        self.mirnas = []
        self.database = database
        self.species_code = species_code
        self.TEMPLATE_SUBSTITUTE_TABLE = {
            "RESOURCE_DIR": resource_dir,
            "SPECIES_CODE": self.species_code,
            "MIRBASE_VERSION": self.database.lower()
            .replace("mirna_", "")
            .replace("mirbase_", "")
            if self.database.startswith("mirna_")
            else self.database,
        }

        self.stemloop_start = self.config["Feature Annotation"]["mirna_stemloop_start"]
        self.stemloop_end = self.config["Feature Annotation"]["mirna_stemloop_end"]

        if self.config["Database Resource"].get("mirbase_gff"):
            self.mirbase_gff = Template(
                self.config["Database Resource"]["mirbase_gff"]
            ).safe_substitute(self.TEMPLATE_SUBSTITUTE_TABLE)

        host = host if host else self.config["Database Resource"].get("mirbase_host")
        user = user if user else self.config["Database Resource"].get("mirbase_user")
        password = (
            password
            if password
            else self.config["Database Resource"].get("mirbase_password")
        )
        if host:
            try:
                self.db_connection = sql.connect(
                    host=host,
                    database=self.database,
                    user=user,
                    password=password,
                )
                self.db_cursor = self.db_connection.cursor()
            except Exception as e:
                errmsg = "Failed to connect to database {} at server {} : {}".format(
                    self.database, host, e
                )
                logging.critical(errmsg)
                raise Exception(errmsg)

        if all(not hasattr(self, x) for x in ["mirbase_gff", "db_connection"]):
            errmsg = "miRBase Connection is not configured properly!"
            logging.critical(errmsg)
            raise Exception(errmsg)

    def get_species_codes(self):
        species_codes = None
        if hasattr(self, "db_cursor"):
            self.db_cursor.execute(self.config["Canned Query"]["mirna_species_query"])
            species_codes = [x[0] for x in self.db_cursor.fetchall()]

        return species_codes

    def get_mirnas(self):

        if not self.mirnas and hasattr(self, "db_connection"):
            self.mirnas = self.load_mirnas_from_db()
        elif not self.mirnas and hasattr(self, "mirbase_gff"):
            self.mirnas = self.load_mirnas_from_gff(self.mirbase_gff)

        return self.mirnas

    def get_mirnas_by_id(self, mirna_id):
        mirnas = []
        for mirna_t in self.get_mirnas():
            if mirna_t.get("mirna_id") == mirna_id:
                mirnas.append(mirna_t)

        return mirnas

    def get_mature_mirnas_by_id(self, mature_acc, mirna_id=None):

        if mirna_id:
            mirnas = []
            mirna_t = self.get_mirnas_by_id(mirna_id)
            if mirna_t:
                mirnas.append(mirna_t)
        else:
            mirnas = self.get_mirnas()

        mature_mirnas = []
        for mirna_t in mirnas:
            if not mirna_t.get("mature_mirnas"):
                continue
            for x in mirna_t.get("mature_mirnas"):
                if x.get("mature_acc") == mature_acc:
                    mature_mirnas.append(x)

        return mature_mirnas

    def get_mirna_ids(self):

        mirna_ids = list(set([x.get("mirna_id") for x in self.get_mirnas()]))

        return mirna_ids

    def get_mature_mirna_acc(self):
        mature_mirna_accs = (x.get("mature_acc") for x in self.get_mirna_mimat_ids())

        return list(set(mature_mirna_accs))

    def get_mirna_mimat_ids(self):

        mirna_mimat_ids = []
        for mirna in self.get_mirnas():
            if mirna.get("mature_mirnas"):
                for x in mirna.get("mature_mirnas"):
                    mirna_mimat_ids.append(
                        {
                            "mirna_id": mirna.get("mirna_id"),
                            "mature_acc": x.get("mature_acc"),
                        }
                    )

        return mirna_mimat_ids

    def load_mirnas_from_db(self):
        """
        Get mirnas from miRBase MySQL database

        """
        logging.info("Loading miRNAs from database {}".format(self.database))
        mirna_query = Template(
            self.config["Canned Query"]["mirna_query"]
        ).safe_substitute(self.TEMPLATE_SUBSTITUTE_TABLE)
        self.db_cursor.execute(mirna_query)
        column_names = [i[0] for i in self.db_cursor.description]
        mirnas = [dict(zip(column_names, x)) for x in self.db_cursor.fetchall()]

        df = (
            pandas.DataFrame(mirnas)
            .groupby(Mirbase.mirna_attrs)
            .agg(list)
            .reset_index()
        )

        self.mirnas = [self.transform_mirna_db_rec(x) for x in df.to_dict("records")]
        return self.mirnas

    @staticmethod
    def transform_mirna_db_rec(mirna):
        mature_mirnas = [
            dict(zip(Mirbase.mature_mirna_attrs, y))
            for y in zip(*[mirna.get(x) for x in Mirbase.mature_mirna_attrs])
        ]

        for mature_mirna in mature_mirnas:
            mature_from = int(mature_mirna.get("mature_from"))
            mature_to = int(mature_mirna.get("mature_to"))
            if mirna.get("strand") == "+":
                mature_mirna["start"] = mirna.get("start") + mature_from - 1
                mature_mirna["end"] = mirna.get("start") + mature_to - 1
            else:
                mature_mirna["start"] = mirna.get("end") - mature_to + 1
                mature_mirna["end"] = mirna.get("end") - mature_from + 1

        mirna["mature_mirnas"] = mature_mirnas
        for field in Mirbase.mature_mirna_attrs:
            if field in mirna:
                del mirna[field]
        return mirna

    def load_mirnas_from_gff(self, gff=None):
        """

        Args:
            gff: GFF file from miRBase website

        Returns:

        """

        required_gff_fields = ["Derives_from", "Alias", "ID"]

        if not gff:
            if hasattr(self, "mirbase_gff"):
                gff = getattr(self, "mirbase_gff")
            else:
                logging.error("Missing miRBase GFF to load")
                return self.mirnas

        logging.info("Loading miRNAs from GFF: {}".format(gff))
        records = gff_to_records(gff)
        if not records:
            errmsg = "Failure loading Mirbase GFF {}".format(gff)
            logging.critical(errmsg)
            raise Exception(errmsg)
        df = (
            pandas.DataFrame(records)
            .astype({"start": "int64", "end": "int64"})
            .rename(
                columns={"derives_from": "Derives_from", "accession_number": "Alias"}
            )
        )
        if any(x not in df.columns for x in required_gff_fields):
            errmsg = "Mirbase GFF {} is missing required fields {}.".format(
                gff, required_gff_fields
            )
            logging.critical(errmsg)
            raise Exception(errmsg)

        df["Derives_from"] = df.apply(
            lambda x: df.loc[
                (
                    (df["Alias"] == x["Derives_from"])
                    | (df["ID"] == x["Derives_from"])
                    | (df["type"] == "mirna_primary_transcript")
                )
                & (df["start"] <= x["start"])
                & (df["end"] >= x["end"])
            ].iloc[0]["ID"]
            if pandas.notna(x["Derives_from"])
            else "",
            axis=1,
        )
        df_primary = df.loc[df["type"] == "miRNA_primary_transcript"]
        df_primary.set_index("ID", drop=False, inplace=True)
        gff_groupby_columns = GFF_COLUMNS[0:8] + ["Derives_from"]
        for item in ["start", "end"]:
            gff_groupby_columns.remove(item)

        df_mature = (
            df.loc[df["type"] == "miRNA"]
            .groupby(gff_groupby_columns)
            .agg(list)
            .reset_index()
        )

        mirbase_gff_attrs_rename = {
            "Name": "mirna_id",
            "Name_mature": "mature_name",
            "Alias_mature": "mature_acc",
            "start_mature": "mature_from",
            "end_mature": "mature_to",
        }

        df_mirnas = df_mature.join(
            df_primary, on="Derives_from", lsuffix="_mature"
        ).rename(columns=mirbase_gff_attrs_rename)[
            Mirbase.mirna_attrs + Mirbase.mature_mirna_attrs
        ]
        self.mirnas = [
            self.transform_mirbase_gff_rec(x) for x in df_mirnas.to_dict("records")
        ]
        return self.mirnas

    @staticmethod
    def transform_mirbase_gff_rec(mirna):
        mature_mirnas = [
            dict(zip(Mirbase.mature_mirna_attrs, y))
            for y in zip(*[mirna.get(x) for x in Mirbase.mature_mirna_attrs])
        ]

        for mature_mirna in mature_mirnas:
            mature_mirna["start"] = int(mature_mirna.get("mature_from"))
            mature_mirna["end"] = int(mature_mirna.get("mature_to"))

        mirna["mature_mirnas"] = mature_mirnas
        for field in Mirbase.mature_mirna_attrs:
            if field in mirna:
                del mirna[field]
        return mirna

    def get_features(self):
        """
        Overrides get_features method of FeatureDatabase

        Get annotated miRNA features

        Returns:
            (list): list of features

        """

        if not self.features:
            features = []
            for mirna in self.get_mirnas():
                features.extend(
                    self.annotate_mirna_feature(
                        mirna,
                        stemloop_start=self.stemloop_start,
                        stemloop_end=self.stemloop_end,
                    )
                )
            self.features = features
        self.update_display_label(self.features)
        return self.features

    @staticmethod
    def annotate_mirna_feature(mirna, stemloop_start=1, stemloop_end=6):
        """
        Determine the precursor, mature and stemloop regions of miRNA:
        The smaller of the two regions flanking the mature mirna is taken as precursor.
        Annotated types: precursor, mature or star, stemloop (if not star)
        The regions that fall outside precursor, mature and star strands and stemloops are "unannotated"

        Args:
            mirna (dict):

        Returns:
            features (list): list of annotated features
        """

        features = []

        mature_mirnas = mirna.get("mature_mirnas")
        for mature_mirna in mature_mirnas:
            if mature_mirna.get("start") - mirna.get("start") < mirna.get(
                "end"
            ) - mature_mirna.get("end"):
                if mirna.get("start") < mature_mirna.get("start"):
                    mature_mirna["pre_start"] = mirna.get("start")
                    mature_mirna["pre_end"] = mature_mirna.get("start") - 1
                if "*" not in mature_mirna.get("mature_name"):
                    mature_mirna["stem_start"] = (
                        mature_mirna.get("end") + stemloop_start
                    )
                    mature_mirna["stem_end"] = mature_mirna.get("end") + stemloop_end
            else:
                if mature_mirna.get("end") < mirna.get("end"):
                    mature_mirna["pre_start"] = mature_mirna.get("end") + 1
                    mature_mirna["pre_end"] = mirna.get("end")
                if "*" not in mature_mirna.get("mature_name"):
                    mature_mirna["stem_start"] = mature_mirna["start"] - stemloop_end
                    mature_mirna["stem_end"] = mature_mirna["start"] - stemloop_start

            features.append(
                {
                    "feature_type": "star"
                    if "*" in mature_mirna.get("mature_name")
                    else "mature",
                    "start": mature_mirna["start"],
                    "end": mature_mirna["end"],
                    "name": "{},{}".format(
                        mirna.get("mirna_id"), mature_mirna.get("mature_acc")
                    ),
                    "chr": mirna.get("chr"),
                    "strand": mirna.get("strand"),
                }
            )
            if mature_mirna.get("pre_start"):
                features.append(
                    {
                        "feature_type": "precursor",
                        "start": mature_mirna["pre_start"],
                        "end": mature_mirna["pre_end"],
                        "name": mirna.get("mirna_id"),
                        "chr": mirna.get("chr"),
                        "strand": mirna.get("strand"),
                    }
                )
            if mature_mirna.get("stem_start"):
                features.append(
                    {
                        "feature_type": "stemloop",
                        "start": mature_mirna["stem_start"],
                        "end": mature_mirna["stem_end"],
                        "name": mirna.get("mirna_id"),
                        "chr": mirna.get("chr"),
                        "strand": mirna.get("strand"),
                    }
                )

        mature_mirnas = sorted(mature_mirnas, key=lambda x: x.get("start"))
        mstart = min(
            [
                mature_mirnas[0].get(x)
                for x in ["stem_start", "pre_start", "start"]
                if mature_mirnas[0].get(x)
            ]
        )
        mend = max(
            [
                mature_mirnas[0].get(x)
                for x in ["stem_end", "pre_end", "end"]
                if mature_mirnas[0].get(x)
            ]
        )

        if len(mature_mirnas) == 1:
            if mstart == mirna.get("start"):
                anon_start = mend + 1
                anon_end = mirna.get("end")
            else:
                anon_start = mirna.get("start")
                anon_end = mstart - 1
        else:
            mstart2 = min(
                [
                    mature_mirnas[1].get(x)
                    for x in ["stem_start", "pre_start", "start"]
                    if mature_mirnas[1].get(x)
                ]
            )
            mend2 = max(
                [
                    mature_mirnas[1].get(x)
                    for x in ["stem_end", "pre_end", "end"]
                    if mature_mirnas[1].get(x)
                ]
            )
            anon_start = mend + 1
            anon_end = mstart2 - 1

        if anon_start <= anon_end:
            features.append(
                {
                    "feature_type": "unannotated",
                    "start": anon_start,
                    "end": anon_end,
                    "name": mirna.get("mirna_id"),
                    "chr": mirna.get("chr"),
                    "strand": mirna.get("strand"),
                }
            )

        return sorted(features, key=lambda x: x.get("start"))

    def save_features(self):
        if hasattr(self, "feature_csv"):
            self.dump_features(self.feature_csv)

    def close(self):
        if hasattr(self, "db_connection"):
            self.db_cursor.close()
            self.db_connection.close()


def get_common_id(mlist):

    prefix = os.path.commonprefix(mlist).rstrip("-")
    return prefix


def main():
    print(timeit.timeit(lambda: annotate("mirna_22"), number=1))


def annotate(database):

    mirbase = Mirbase(database=database, species_code="hsa")
    mirbase.dump_features("mirbase_{}_features.csv".format(database))

    json.dumps(
        mirbase.get_overlapping_features_from_hierarchy("chr1", 17369, 17408, "-"),
        indent=True,
    )

    mirbase.close()


if __name__ == "__main__":
    main()
