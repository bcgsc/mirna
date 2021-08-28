import json
import os

from mirna_profiling.library_stat.alignment_stat import get_project_stat, AlignmentStat


def test_update_nomenclature(mirbase_21_db):

    isofile = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "mock_input", "isoforms.json"
    )
    isomir_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "mock_input",
        "isoforms_w_isomir.json",
    )
    with open(isofile, "r") as jf:
        isoforms = json.load(jf)

    with open(isomir_file, "r") as jf:
        isomirs = json.load(jf)

    isoforms = AlignmentStat.update_isoform_nomenclature(isoforms, mirbase_21_db)

    assert isomirs == isoforms


def test_get_project_stat(hg38_lib_data):
    MIRNA_FILES = ["isoforms.txt", "mirna_species.txt", "miRNA.txt"]
    ucsc = "hg38"
    mirbase_version = "mirna_21"
    species_code = "hsa"

    libname = hg38_lib_data.get("library_name")
    proj_dir = hg38_lib_data.get("library_directory")
    out_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "mock_output")
    proj_stats = get_project_stat(
        proj_dir,
        ucsc,
        mirbase_version,
        species_code,
        out_dir=out_dir,
    )

    assert len(proj_stats) == 1
    lib_stats = [x for x in proj_stats if x.get("Library") == libname]

    bam_stat_file = hg38_lib_data.get("bam_stat_json")
    with open(bam_stat_file, "r") as jf:
        bam_stat_json = json.load(jf)

    assert lib_stats and bam_stat_json.items() <= lib_stats[0].items()

    assert all(
        os.path.exists(os.path.join(out_dir, libname, "{}_features".format(libname), f))
        for f in MIRNA_FILES
    )
