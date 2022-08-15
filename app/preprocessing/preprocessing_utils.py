import json
import os
import pandas as pd
from zipfile import ZipFile


def get_series_probesets(accessions):
    '''Returns a dict of accessions with their respective probesets'''

    with open("json/series_platforms.json") as series_platforms_json:
        series_platforms = json.load(series_platforms_json)

    with open("json/platform_probesets.json") as platform_probesets_json:
        platform_probesets = json.load(platform_probesets_json)

    supported_probesets = set()
    for probemap_fname in os.listdir("probe_maps"):
        fname_base = probemap_fname.split(".")[0]
        probeset = fname_base.split("_", 1)[1]
        supported_probesets.add(probeset)

    series_probesets = dict()

    for accession in accessions:
        accession = accession.lower()
        if accession in series_platforms.keys():
            platform = series_platforms[accession]
        else:
            continue
        if (platform in platform_probesets.keys() and
                platform_probesets[platform] in supported_probesets):
            series_probesets[accession] = platform_probesets[platform]
        for probeset in platform_probesets.values():
            if platform in probeset and probeset in supported_probesets:
                series_probesets[accession] = probeset

    return series_probesets


def get_valid_geo_accessions(accessions):
    '''Returns a list of valid GEO accessions from given list of accessions'''

    valid_accessions = list(get_series_probesets(accessions).keys())

    return valid_accessions


def get_valid_gdc_projects(project_names):
    '''Returns a list of valid GDC projects from project_names'''

    with open("json/gdc_projects.json") as gdc_projects_json:
        gdc_projects = json.load(gdc_projects_json)

    valid_projects = []
    for project_name in project_names:
        if "-" not in project_name:
            continue
        proj_name_lower = project_name.lower()
        proj_base, proj_ext = proj_name_lower.split("-", 1)
        if proj_base in gdc_projects and proj_ext in gdc_projects[proj_base]:
            valid_projects.append(project_name)

    return valid_projects


def combine_gdc(data_dir):
    counts_paths = get_counts_paths(data_dir)
    coldata_paths = get_coldata_paths(data_dir)
    counts_dfs = []
    coldata_dfs = []

    for counts_path, coldata_path in zip(counts_paths, coldata_paths):
        counts_df = pd.read_csv(counts_path, sep="\t")
        coldata_df = pd.read_csv(coldata_path, sep="\t")
        counts_df.set_index("Hugo_Symbol", inplace=True)
        coldata_df.set_index("sample_name", inplace=True)
        counts_dfs.append(counts_df)
        coldata_dfs.append(coldata_df)

    combined_counts = pd.concat(counts_dfs,axis=1)
    combined_coldata = pd.concat(coldata_dfs,axis=0)

    combined_counts.to_csv(os.path.join(data_dir, "counts_processed.tsv"), sep="\t")
    combined_coldata.to_csv(os.path.join(data_dir, "coldata_processed.tsv"), sep="\t")


def get_counts_paths(data_dir):
    counts_paths = []

    for fname in os.listdir(data_dir):
        if "_counts.tsv" in fname:
            counts_paths.append(os.path.join(data_dir, fname))

    return counts_paths


def get_unmapped_counts_paths(data_dir):
    counts_paths = []

    for fname in os.listdir(data_dir):
        if "_counts_unmapped.tsv" in fname:
            counts_paths.append(os.path.join(data_dir, fname))
    
    return counts_paths


def get_coldata_paths(data_dir):
    coldata_paths = []

    for fname in os.listdir(data_dir):
        if "_coldata.tsv" in fname:
            coldata_paths.append(os.path.join(data_dir, fname))
    
    return coldata_paths


def map_probes(counts_path, species):

    counts_fname = counts_path.split("/")[-1]
    accession_id = counts_fname.split("_")[0].lower()
    counts_df = pd.read_csv(counts_path, sep="\t")

    probeset = get_series_probesets([accession_id])[accession_id]
    if not probeset:
        return ""
    counts_df.rename(columns={"probe": probeset}, inplace=True)
    probe_map = pd.read_csv(f"probe_maps/{species}_{probeset}.tsv", sep="\t")
    try:
        counts_df = pd.merge(counts_df, probe_map, on=probeset, how="outer")
    except:
        print("Something went wrong while mapping probes to gene symbols.")

    cols = list(counts_df.columns)
    cols = [cols[-1]] + cols[1:-1]
    counts_df = counts_df[cols]
    counts_df.dropna(inplace=True)

    return counts_df


def prep_geo_counts(data_dir):
    '''Prepares unmapped expression matrices from GEO'''

    counts_paths = get_unmapped_counts_paths(data_dir)
    for counts_path in counts_paths:
        counts_df = map_probes(counts_path, "hsapiens")
        new_counts_path = counts_path.replace("_unmapped.tsv", "_processed.tsv")
        counts_df.to_csv(new_counts_path, sep="\t", index=False)

    # Cleanup unmapped counts
    for counts_path in counts_paths:
        os.remove(counts_path)


def zip_preprocessed_data(data_dir):
    '''Zips counts and coldata files and returns filename'''

    old_wd = os.getcwd()
    new_wd = data_dir
    os.chdir(new_wd)

    files_to_zip = []
    is_geo = False
    for fname in os.listdir():
        if fname.endswith("processed.tsv"):
            files_to_zip.append(fname)
            is_geo = True if "gse" in fname.lower() else False

    source = "geo" if is_geo else "gdc"
    zip_fname = f"{source}_processed.zip"
    print(zip_fname)
    with ZipFile(zip_fname, "w") as zip_file:
        for file in files_to_zip:
            print(f"zipping file: {file}")
            zip_file.write(file)

    # clean up files
    for fpath in files_to_zip:
        os.remove(fpath)

    os.chdir(old_wd)

    return zip_fname
