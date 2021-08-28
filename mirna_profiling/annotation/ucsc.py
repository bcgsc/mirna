import json
import logging
import os
import re
import timeit
from string import Template

# Password authentication problem with mysql-connector-python for MySQL 5.5 and older
# Version mismatch ("Bad handshake") with mysql-connector-python for MySQL 5.0
# import mysql.connector as sql
import pandas
import pymysql as sql

from mirna_profiling.annotation.feature_db import FeatureDatabase


class UCSC(FeatureDatabase):
    """
    A class that acts as an interface to UCSC, and retrieves features of interest from the UCSC MySQL databases

    """

    def __init__(
        self,
        database,
        host=None,
        user=None,
        password=None,
        **kwarg,
    ):
        """
        Initialize database connection

        Raises:
            Exception if database connection fails

        """
        super().__init__(**kwarg)
        # connect_string = 'mysql://genome@genome-mysql.soe.ucsc.edu:3606/hg38'
        self.database = database
        host = host if host else self.config["Database Resource"].get("ucsc_host")
        user = user if user else self.config["Database Resource"].get("ucsc_user")
        password = (
            password
            if password
            else self.config["Database Resource"].get("ucsc_password")
        )

        if host:
            try:
                self.db_connection = sql.connect(
                    host=host,
                    database=database,
                    user=user,
                    password=password,
                )
                self.db_cursor = self.db_connection.cursor()
            except Exception as e:
                errmsg = (
                    "Failed to connect to UCSC database {} at server {}: {}".format(
                        self.database, host, e
                    )
                )
                logging.critical(errmsg)
                raise Exception(errmsg)

        if not hasattr(self, "db_connection"):
            errmsg = "UCSC Database is not configured properly!"
            logging.critical(errmsg)
            raise Exception(errmsg)

        self.repeats = []
        self.genes = []

    def get_db_tables(self, pattern=None):
        """
        Get tables whose name match a pattern

        Args:
            pattern (str): pattern of table name to query

        Returns:
            tables (list): list of table names matching the pattern specified
        """
        query = "SHOW TABLES"
        if pattern:
            query = "{} like {}".format(query, pattern)
        self.db_cursor.execute(query)
        tables = [x[0] for x in self.db_cursor.fetchall()]
        return tables

    def load_repeats_from_db(self):
        """
        Get repeat elements from UCSC RepeatMasker tables.

        For most genomes, there is only one table, "rmsk", for all chromosomes,
        but for some genomes, e.g., mm9, there is one table per chromosome, e.g. chr1_rmsk, or chr1_random_rmsk

        Returns:
            (list of dict): repeats on canonical chrs (i.e, those with one _ in the table name)

        """
        repeat_tables = [x for x in self.get_db_tables("'%rmsk'") if x.count("_") <= 1]
        repeats = []
        logging.info("Loading UCSC repeats from database {}".format(self.database))
        for repeat_table in repeat_tables:
            repeat_query = Template(
                self.config["Canned Query"]["repeat_query"]
            ).safe_substitute({"REPEAT_TABLE": repeat_table})
            self.db_cursor.execute(repeat_query)
            column_names = [i[0] for i in self.db_cursor.description]
            for row in self.db_cursor.fetchall():
                rec = dict(zip(column_names, row))
                if all(
                    not re.match(y, rec.get("genoName"))
                    for y in self.config["Canned Query"]["ucsc_chr_exclude_pattern"]
                ):
                    repeats.append(rec)
        return repeats

    def load_genes_from_db(self, gene_table=None):
        """
        Get (Known or Refseq) genes from gene table

        Args:
            gene_table (str): name of gene table to load. Default: knownGene or refGene

        Returns:
            (list of dict): list of genes
        """
        if not gene_table:
            gtables = self.get_db_tables("'knownGene'")
            if not gtables:
                gtables = self.get_db_tables("'refGene'")
            gene_table = gtables[0] if gtables else None

        if not gene_table:
            logging.warning("Missing gene_table for query")
            return []

        logging.info(
            "Loading UCSC genes from database {} table {}".format(
                self.database, gene_table
            )
        )
        gene_query = Template(
            self.config["Canned Query"]["gene_query"]
        ).safe_substitute(
            {
                "GENE_TABLE": gene_table,
                "GENE_FK": "refseq" if gene_table == "refGene" else "kgID",
            }
        )
        genes = [self.format_gene(x) for x in self.execute_query(gene_query)]

        return genes

    def get_genes(self):
        if not self.genes and hasattr(self, "db_connection"):
            self.genes = self.load_genes_from_db()
        return self.genes

    def get_repeats(self):
        if not self.repeats and hasattr(self, "db_connection"):
            self.repeats = self.load_repeats_from_db()

        return self.repeats

    def get_features(self):
        """
        Overrides get_features method of FeatureDatabase

        Get gene and repeat features

        Returns:
            (list): list of features

        """
        if not self.features:
            if hasattr(self, "feature_csv") and os.path.exists(
                getattr(self, "feature_csv")
            ):
                self.features = self.load_features_from_csv(
                    getattr(self, "feature_csv")
                )
            else:
                self.features = []
                for gene in self.get_genes():
                    self.features.extend(self.annotate_gene_features(gene))
                self.features.extend(
                    [self.annotate_repeat_feature(x) for x in self.get_repeats()]
                )
        self.update_display_label(self.features)
        return self.features

    def load_kgxref_from_db(self):
        """
        Get gene xrefs, with geneSymbol, Refseq and UCSC gene id mapping

        Returns:

        """
        return self.execute_query(self.config["Canned Query"]["XREF_QUERY"])

    def execute_query(self, query):
        """
        Execute SQL query and returns results

        Args:
            query (str):

        Returns:
            (list of dict): query results
        """
        results = []
        self.db_cursor.execute(query)
        column_names = [i[0] for i in self.db_cursor.description]
        results = [dict(zip(column_names, x)) for x in self.db_cursor.fetchall()]
        return results

    def close(self):
        if hasattr(self, "db_connection"):
            self.db_cursor.close()
            self.db_connection.close()

    def format_gene(self, gene):
        """

        Args:
            gene (dict):

        Returns:
            gene (dict): formatted gene
        """
        for key in gene.keys():
            if isinstance(gene.get(key), bytes):
                gene[key] = gene[key].decode("UTF-8")

        # add 1 to 0 based ucsc coordinates
        # ucsc database uses 0-based start and 1-based end coordinates
        # so add 1 to start only, but not end
        # http://genome.ucsc.edu/FAQ/FAQtracks#tracks1
        gene["txStart"] += 1
        if gene["cdsStart"] != gene["cdsEnd"]:
            gene["cdsStart"] += 1
        gene["exonStarts"] = ",".join(
            [str(int(x) + 1) if x else x for x in gene["exonStarts"].split(",")]
        )
        return gene

    @staticmethod
    def annotate_repeat_feature(repeat):
        feature = {
            "chr": repeat.get("genoName"),
            "strand": repeat.get("strand"),
            "start": repeat.get("genoStart") + 1,
            "end": repeat.get("genoEnd"),
            "feature_type": repeat.get("repClass"),
            "name": repeat.get("repName"),
        }

        return feature

    @staticmethod
    def annotate_gene_features(gene):
        """
        Annotate features in UCSC gene region
        Args:
            gene (dict):

        Returns:
            features (list):

        """
        features = []

        exonstarts = [int(x) for x in gene.get("exonStarts").rstrip(",").split(",")]
        exonends = [int(x) for x in gene.get("exonEnds").rstrip(",").split(",")]

        exons = zip(exonstarts, exonends)

        if gene.get("cdsStart") == gene.get("cdsEnd"):
            if gene.get("geneSymbol").startswith(
                "SNOR"
            ) or "small nucleolar RNA" in gene.get("description"):
                feature_type = "snoRNA"
            else:
                feature_type = "Non-Coding Exon"

            for (exonstart, exonend) in exons:
                features.append(
                    {
                        "name": gene.get("geneSymbol"),
                        "feature_type": feature_type,
                        "start": exonstart,
                        "end": exonend,
                        "chr": gene.get("chrom"),
                        "strand": gene.get("strand"),
                    }
                )
        else:

            for (exonstart, exonend) in exons:
                utr5_start = utr5_end = utr3_start = utr3_end = None
                coding_start = coding_end = None
                if exonstart < gene.get("cdsStart"):
                    utr5_start = exonstart
                    utr5_end = (
                        gene.get("cdsStart") - 1
                        if exonend >= gene.get("cdsStart")
                        else exonend
                    )
                if exonend > gene.get("cdsEnd"):
                    utr3_end = exonend
                    utr3_start = (
                        gene.get("cdsEnd") + 1
                        if exonstart <= gene.get("cdsEnd")
                        else exonstart
                    )

                if exonend >= gene.get("cdsStart") and exonstart <= gene.get("cdsEnd"):
                    coding_start = (
                        exonstart
                        if exonstart > gene.get("cdsStart")
                        else gene.get("cdsStart")
                    )
                    coding_end = (
                        exonend if exonend < gene.get("cdsEnd") else gene.get("cdsEnd")
                    )

                if utr5_start and utr5_end:
                    features.append(
                        {
                            "name": gene.get("geneSymbol"),
                            "feature_type": "5 UTR"
                            if gene.get("strand") == "+"
                            else "3 UTR",
                            "start": utr5_start,
                            "end": utr5_end,
                            "chr": gene.get("chrom"),
                            "strand": gene.get("strand"),
                        }
                    )
                if coding_start and coding_end:
                    features.append(
                        {
                            "name": gene.get("geneSymbol"),
                            "feature_type": "Coding Exon",
                            "start": coding_start,
                            "end": coding_end,
                            "chr": gene.get("chrom"),
                            "strand": gene.get("strand"),
                        }
                    )
                if utr3_start and utr3_end:
                    features.append(
                        {
                            "name": gene.get("geneSymbol"),
                            "feature_type": "3 UTR"
                            if gene.get("strand") == "+"
                            else "5 UTR",
                            "start": utr3_start,
                            "end": utr3_end,
                            "chr": gene.get("chrom"),
                            "strand": gene.get("strand"),
                        }
                    )

        for idx in range(0, len(exonstarts) - 1):
            features.append(
                {
                    "name": gene.get("geneSymbol"),
                    "feature_type": "Intron",
                    "start": exonends[idx] + 1,
                    "end": exonstarts[idx + 1] - 1,
                    "chr": gene.get("chrom"),
                    "strand": gene.get("strand"),
                }
            )

        return sorted(features, key=lambda x: x.get("start"))


def main():
    print(timeit.timeit(lambda: annotate("hg38"), number=1))


def annotate(genome):

    ucsc = UCSC(database=genome)
    genes = ucsc.get_genes()
    pandas.DataFrame(genes).to_csv(
        "ucsc_{}_genes.csv".format(genome), sep="\t", index=False
    )
    ucsc.dump_features("ucsc_{}_features.csv".format(genome))

    json.dumps(
        ucsc.get_overlapping_features_from_hierarchy("chr1", 13221, 14409, "+"),
        indent=True,
    )

    ucsc.close()


if __name__ == "__main__":
    main()
