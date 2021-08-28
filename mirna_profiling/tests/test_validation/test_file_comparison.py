import os
from mirna_profiling.validation.file_comparison import compare_directory, compare_csv


def test_compare_directory(hg38_lib_data):
    libname = hg38_lib_data.get("library_name")
    feature_dir = os.path.join(
        hg38_lib_data.get("library_directory"), "{}_features".format(libname)
    )
    dir2 = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "mock_input", libname
    )
    assert not compare_directory(feature_dir, dir2)


def test_compare_csv(hg38_lib_data):

    libname = hg38_lib_data.get("library_name")
    feature_dir = os.path.join(
        hg38_lib_data.get("library_directory"), "{}_features".format(libname)
    )

    dir2 = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "mock_input", libname
    )

    file1 = os.path.join(feature_dir, "filtered_taglengths.csv")
    file2 = os.path.join(dir2, "filtered_taglengths.csv")
    assert compare_csv(file1, file2)

    file1 = os.path.join(feature_dir, "softclip_taglengths.csv")
    file2 = os.path.join(dir2, "softclip_taglengths.csv")
    assert not compare_csv(file1, file2)
