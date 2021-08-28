import pytest
import os

from mirna_profiling.annotation.mirbase import Mirbase

mock_input_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "mock_input")
mock_output = os.path.join(os.path.dirname(os.path.realpath(__file__)), "mock_output")
lib_hg38 = "m00038"


@pytest.fixture(scope="module")
def mirbase_21_db():

    mirbase = Mirbase("mirna_21", "hsa")
    return mirbase


@pytest.fixture(scope="module")
def hg38_proj_dir():

    return os.path.join(mock_input_dir, "hg38")


@pytest.fixture(scope="module")
def hg38_lib_data(hg38_proj_dir):
    hg38_lib_dir = os.path.join(hg38_proj_dir, lib_hg38)
    lib_data = {
        "library_name": lib_hg38,
        "library_directory": hg38_lib_dir,
        "bam": os.path.join(hg38_lib_dir, "{}.bam".format(lib_hg38)),
        "annotated_bam": os.path.join(hg38_lib_dir, "{}.bam.annot".format(lib_hg38)),
        "adapter_report": os.path.join(
            hg38_lib_dir, "{}_adapter.report".format(lib_hg38)
        ),
        "bam_stat_json": os.path.join(
            hg38_lib_dir, "{}_bam_stat.json".format(lib_hg38)
        ),
    }

    return lib_data
