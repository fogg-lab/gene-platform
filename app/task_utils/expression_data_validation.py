from typing import List


def check_coldata_has_required_columns(coldata_rows: List[List[str]], require_batch=False) -> str:
    """
    Generates an error message if coldata doesn't have a required column.
    Args:
        coldata_rows List[List[str]]: Contents of coldata.tsv.
    Returns:
        str: Error message if coldata does not have a required column.
    """

    required_columns = ["sample_name", "condition"]
    if require_batch:
        required_columns.append("batch")

    for required_column in required_columns:
        if required_column not in coldata_rows[0]:
            return f"coldata file does not have a '{required_column}' column"

    return ""


def check_counts_matches_coldata(counts_colnames: List[str], coldata_rows: List[List[str]]) -> str:
    """
    Generates an error message if counts doesn't have all samples listed in coldata.
    Args:
        counts_colnames (List[str]): List of column names from counts file.
        coldata_rows (List[List[str]]): Contents of coldata.tsv.
    Returns:
        str: Error message if sample names in coldata and counts files don't match.
    """

    coldata_sample_name_idx = coldata_rows[0].index("sample_name")
    coldata_sample_names = [row[coldata_sample_name_idx] for row in coldata_rows[1:]]

    samples = []

    coldata_rows.pop(0)
    for coldata_row in coldata_rows:
        if coldata_row:
            samples.append(coldata_row[0])

    for sample in coldata_sample_names:
        if sample not in counts_colnames:
            return f"Sample '{sample}' in coldata not found in counts"


def check_factor_levels(reference_level: str, contrast_level: str,
                        coldata_rows: List[List[str]]) -> str:
    """
    Generates an error message if factor levels aren't present in the coldata file.
    Args:
        reference_level (str): Reference level specified in the config parameters.
        contrast_level (str): Contrast level specified in the config parameters.
        coldata_rows (List[List[str]]): Contents of coldata.tsv.
    Returns:
        str: Error message if factor levels aren't present in the coldata.
    """

    condition_col_index = 0
    col_header_row = [colname.lower() for colname in coldata_rows[0]]
    if "condition" in col_header_row:
        condition_col_index = col_header_row.index("condition")
    else:
        return "'condition' column not present in coldata (line 1)\n"

    contrast_level_found = False
    reference_level_found = False

    # remove header row from coldata
    coldata_rows.pop(0)
    err_msg = ""
    for coldata_row in coldata_rows:
        if coldata_row:
            factor_level = coldata_row[condition_col_index]
            if factor_level == contrast_level:
                contrast_level_found = True
            elif factor_level == reference_level:
                reference_level_found = True
            else:
                return f"Unknown factor level '{factor_level}'"
    if contrast_level_found is False:
        err_msg += f"Unknown contrast level '{contrast_level}'\n"
    if reference_level_found is False:
        err_msg += f"Unknown reference level '{reference_level}'\n"

    return err_msg
