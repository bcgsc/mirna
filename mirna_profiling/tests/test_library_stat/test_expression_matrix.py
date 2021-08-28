import os

from mirna_profiling.library_stat.expression_matrix import (
    get_mirna_expression_matrix,
    get_mimat_expression_matrix,
    exp_to_csv,
)
from mirna_profiling.validation.file_comparison import compare_files
from mirna_profiling.annotation.mirbase import Mirbase


def test_get_mirna_expression_matrix(mirbase_21_db, hg38_proj_dir):
    EXP_FILES = ["expn_matrix_norm_log.txt", "expn_matrix_mimat_norm_log.txt"]

    out_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        "mock_output",
    )
    os.makedirs(out_dir, 0o775, exist_ok=True)
    exp_df = get_mirna_expression_matrix(hg38_proj_dir, mirbase_21_db)
    exp_mimat_df = get_mimat_expression_matrix(hg38_proj_dir, mirbase_21_db)
    assert exp_df.shape == (1881, 1)
    assert exp_mimat_df.shape == (2588, 1)

    exp_to_csv(exp_df, "expn_matrix", out_dir)
    exp_to_csv(exp_mimat_df, "expn_matrix_mimat", out_dir)

    for exp_f in EXP_FILES:
        file_expected = os.path.join(hg38_proj_dir, exp_f)
        file_actual = os.path.join(out_dir, exp_f)
        assert compare_files(
            file_expected,
            file_actual,
            file_types={"csv": EXP_FILES},
        )
