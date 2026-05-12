import csv
import sys


def set_csv_field_size_limit() -> None:
    """
    Set the maximum field size limit for reading CSV files.

    Notes
    -----
    This function should be called before reading any CSV files to ensure that the field size
    limit is properly set.

    """
    max_int = sys.maxsize

    while max_int > 1:
        try:
            csv.field_size_limit(max_int)
            break
        except OverflowError:
            max_int = int(max_int / 10)
