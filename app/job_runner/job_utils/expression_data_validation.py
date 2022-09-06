def check_factor_levels(reference_level, contrast_level, coldata):
    """
    Ensures factor levels are present in the coldata file
    If factor levels are present, returns empty string
    Otherwise, returns an error message
    Args:
        reference_level: string
        contrast_level: string
        coldata: 2d array
    Returns:
        string
    """

    condition_col_index = 0
    col_header_row = [colname.lower() for colname in coldata[0]]
    if "condition" in col_header_row:
        condition_col_index = col_header_row.index("condition")
    else:
        return "'condition' column not present in coldata (line 1)\n"

    contrast_level_found = False
    reference_level_found = False

    # remove header row from coldata
    coldata.pop(0)
    err_msg = ""
    for coldata_row in coldata:
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


def check_coldata_matches_counts(counts_colnames, coldata_rows):
    """
    Ensure rows in coldata match with the column names for samples in counts
    Assumes that sample names are listed on first row (header) of counts file
    Also assumes that sample names are listed in first column of coldata file,
        starting on second row of coldata file (first row after the header)
    Returns empty string if coldata and counts sample names match
    Otherwise, returns an error message
    Args:
        counts_colnames: array
        coldata_rows: 2d array
    Returns:
        string
    """

    samples = []

    coldata_rows.pop(0)
    for coldata_row in coldata_rows:
        if coldata_row:
            samples.append(coldata_row[0])

    if not samples:
        return "no samples are listed in the coldata file"

    # remove leading elements from counts which are not sample names
    while counts_colnames and counts_colnames[0] != samples[0]:
        counts_colnames.pop(0)

    if not counts_colnames:
        return f"sample '{samples[0]}' from coldata not found in counts file"

    # if counts and coldata match, err_msg will be empty
    status_msg = ""

    # make sure each sample name matches between counts and coldata
    while samples and counts_colnames and samples[0] == counts_colnames[0]:
        samples.pop(0)
        counts_colnames.pop(0)

    if counts_colnames and not samples:
        status_msg = f"sample '{counts_colnames[0]}' not found in coldata file"

    elif samples and not counts_colnames:
        status_msg = f"sample '{samples[0]}' not found in counts file"

    elif samples and counts_colnames and samples[0] != counts_colnames[0]:
        status_msg = f"sample '{samples[0]}' in coldata file does not match the \
            corresponding sample name '{counts_colnames[0]}' in counts file"

    return status_msg
