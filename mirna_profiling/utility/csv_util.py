import pandas
import json
import logging


def csv_to_json(in_csv, out_json=None):
    """
    Convert CSV file to JSON

    Args:
        in_csv:
        out_json:

    Returns:

    """
    try:
        records = pandas.read_csv(in_csv).to_dict("records")
        json_records = json.dumps(records)
        if out_json:
            with open(out_json, "w") as jf:
                json.dump(jf, json_records)
    except Exception as e:
        logging.warning("Error converting CSV {} to JSON".format(in_csv, e))
        json_records = None

    return json_records


def combine_csv(
    csvs,
    axis,
    header=None,
    index_col=None,
    out_csv=None,
    sort_index=False,
    sort_column=False,
    sep=",",
):
    """
    Combine a list of CSV files into one.

    Args:
        out_csv (str): out CSV path
        csvs (list): list of CSV files
        axis: axis of to combine on, 0 for row, 1 for column
        header (int or list of int): row numbers to use as column names
        index_col (str or list of str):

    Returns:

    """

    combined_df = None
    for csv_file in csvs:
        df = pandas.read_csv(csv_file, header=header, sep=sep)
        if index_col:
            df.drop_duplicates(subset=index_col, inplace=True)
            df.set_index(index_col, inplace=True)

        if combined_df is None:
            combined_df = df
        else:
            combined_df = pandas.concat([combined_df, df], axis=axis)
            del df

    if combined_df is not None:
        if sort_index:
            combined_df.sort_index(inplace=True)
        if sort_column:
            combined_df.sort_index(axis=1, inplace=True)

        if out_csv:
            combined_df.to_csv(
                out_csv,
                index=True if index_col is not None else False,
                header=True if header is not None else False,
                sep=sep,
            )

    return combined_df
