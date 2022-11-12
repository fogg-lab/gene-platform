import sys
import os
import subprocess
import argparse
from pathlib import Path
import yaml
from pandas import merge, read_csv, DataFrame

# Import function from helper module located in the parent directory
sys.path.insert(0, str(Path(__file__).parent.parent))
from helper import get_series_probesets


def main():
    args = get_args()

    with open(args.cfg_path, "r") as cfg_file:
        cfg = yaml.safe_load(cfg_file)

    dsets = cfg.get("dsets")
    data_source = cfg.get("data_source")
    if "geo" in data_source.lower():
        subprocess.check_call(f"{args.prep_script} {args.out_dir} "
                          f"{' '.join(dsets)} >>{args.log} 2>&1", shell=True)
    else:
        subprocess.check_call(f"{args.prep_script} {args.data_dir} "
                          f"{args.out_dir} {' '.join(dsets)} >>{args.log} 2>&1", shell=True)
    if "geo" in data_source.lower():
        # Map probes to gene symbols. data_dir for the counts is the output directory.
        prep_geo_counts(args.out_dir, args.probe_maps_dir, args.data_dir, args.out_dir)
    print("Done preprocessing data.")

def get_args():
    """
    Parse command line arguments:
        - Path to data prep script
        - Path to probe map directory (only required if data source is GEO)
        - Path to the persistent data directory for the data source
        - Path to input directory
        - Path to output directory
        - Path to log file
    """
    parser = argparse.ArgumentParser(description="Preprocess data from GEO or GDC.")
    parser.add_argument("-s", "--prep_script", type=str, required=True,
                        help="Path to data prep script, e.g. prep_gdc.r.")
    parser.add_argument("-c", "--cfg_path", type=str, required=True,
                        help="Path to the config file.")
    parser.add_argument("-d", "--data_dir", type=str, required=True,
                        help="Path to persistent data directory for the data source.")
    parser.add_argument("-o", "--out_dir", type=str, required=True,
                        help="Path to the output directory.")
    parser.add_argument("-l", "--log", type=str, required=True,
                        help="Path to log file.")
    parser.add_argument("-p", "--probe_maps_dir", type=str,
                        help="Path to probe map directory (for GEO only).")

    return parser.parse_args()


def prep_geo_counts(data_dir, probe_maps_dir, metadata_dir, out_dir):
    """Prepares unmapped expression matrices from GEO"""

    counts_paths = get_unmapped_counts_paths(data_dir)
    counts_unmapped = set()
    for counts_path in counts_paths:
        counts_df = map_probes(counts_path, "hsapiens", probe_maps_dir, metadata_dir)
        if not isinstance(counts_df, DataFrame) or len(counts_df) == 0:
            counts_unmapped.add(counts_path)
            continue
        counts_fname = Path(counts_path).name
        new_counts_fname = counts_fname.replace("_unmapped.tsv", "_processed.tsv")
        new_counts_path = Path(out_dir) / new_counts_fname
        counts_df.to_csv(new_counts_path, sep="\t", index=False)

    # Cleanup unmapped counts
    for counts_path in counts_paths:
        os.remove(counts_path)

    return counts_unmapped


def map_probes(counts_path, species, probe_map_dir, metadata_dir):
    """Maps probes to gene symbols"""

    counts_fname = Path(counts_path).name
    accession_id = counts_fname.split("_")[0].lower()
    counts_df = read_csv(counts_path, sep="\t")

    probeset = get_series_probesets([accession_id], metadata_dir)[accession_id]
    if not probeset:
        return ""
    counts_df.rename(columns={"probe": probeset}, inplace=True)
    probeset_map_path = Path(probe_map_dir) / f"{species}_{probeset}.tsv"
    if not os.path.exists(probeset_map_path):
        print(f"No probe map found: {probeset_map_path}")
    else:
        print(f"Probe map found: {probeset_map_path}")
    probe_map = read_csv(probeset_map_path, sep="\t",
                            dtype={"hgnc_symbol": object, probeset: object})

    counts_df[probeset] = counts_df[probeset].apply(lambda x: str(x).rstrip("_at"))
    probe_map[probeset] = probe_map[probeset].apply(lambda x: str(x).rstrip("_at"))

    counts_df = merge(counts_df, probe_map, on=probeset, how="outer")

    cols = list(counts_df.columns)
    cols = [cols[-1]] + cols[1:-1]
    counts_df = counts_df[cols]
    counts_df.dropna(inplace=True)

    return counts_df


def get_unmapped_counts_paths(data_dir):
    """Returns a list of paths to unmapped expression matrices"""
    counts_paths = []

    for fname in os.listdir(data_dir):
        if "_counts_unmapped.tsv" in fname:
            counts_paths.append(os.path.join(data_dir, fname))

    return counts_paths


if __name__ == "__main__":
    main()
