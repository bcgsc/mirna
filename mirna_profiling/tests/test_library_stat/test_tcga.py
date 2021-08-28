import os

from mirna_profiling.library_stat.tcga import write_tcga_files
from mirna_profiling.validation.file_comparison import compare_files


def test_write_tcga_files(mirbase_21_db, hg38_lib_data):
    TCGA_FILES = ["mirnas.txt", "isoforms.txt"]
    libname = hg38_lib_data.get("library_name")
    hg38_lib_dir = hg38_lib_data.get("library_directory")

    out_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "mock_output", "tcga"
    )

    write_tcga_files(hg38_lib_dir, "hg38", mirbase_21_db, out_dir=out_dir)

    for test_file in TCGA_FILES:

        file_expected = os.path.join(
            hg38_lib_dir, "{}_features".format(libname), "tcga", test_file
        )
        file_actual = os.path.join(out_dir, libname, test_file)
        assert compare_files(
            file_expected,
            file_actual,
            file_types={"csv": TCGA_FILES},
        )
