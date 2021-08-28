import os

from mirna_profiling.library_stat.adapter import (
    ADAPTER_CATEGORIES,
    ADAPTERS,
    get_bin,
    get_adapter_stat,
)


def test_get_bin():
    assert get_bin(0, ADAPTER_CATEGORIES) == 0
    assert get_bin(14, ADAPTER_CATEGORIES) == 14
    assert get_bin(15, ADAPTER_CATEGORIES) == 25
    assert get_bin(35, ADAPTER_CATEGORIES) == 35
    assert get_bin(1000, ADAPTER_CATEGORIES) == ADAPTER_CATEGORIES[-1]


def test_get_adapter_stat(hg38_lib_data):

    adapter_file = hg38_lib_data.get("adapter_report")
    assert (x in ADAPTERS.values() for x in get_adapter_stat(adapter_file))
