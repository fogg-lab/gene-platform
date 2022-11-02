import json
import os
from pathlib import Path
from zipfile import ZipFile
import pandas as pd

from app.task_utils.task_runner import TaskRunner


class PreprocessingRunner(TaskRunner):
    """Class for preparing and running a preprocessing task."""
    def __init__(self, task_id, task_dir):
        super().__init__(task_id, task_dir)
        self.task_type = "preprocessing"
        self._input_filenames = ["config.yml"]

    def execute_task(self):
        """Run a preprocessing task"""
        return dict(status="", warnings=[], errors=[])

    def validate_config(self, config) -> dict:
        """Ensures config parameters are valid for the task"""
        return dict(status="", warnings=[], errors=[])

    def validate_task(self) -> dict:
        """Validates all input files for the task"""
        return dict(status="", warnings=[], errors=[])

    def _get_unmapped_counts_paths(self):
        """Returns a list of paths to unmapped expression matrices"""
        counts_paths = []

        data_dir = os.path.join(self._task_dir, "data")

        for fname in os.listdir(data_dir):
            if "_counts_unmapped.tsv" in fname:
                counts_paths.append(os.path.join(data_dir, fname))

        return counts_paths

    @staticmethod
    def _get_valid_geo_accessions(accessions):
        """
        Returns a list of valid GEO accessions from given list of accessions
        Args:
            accessions (list[str]): List of GEO accessions
        Returns:
            list: List of valid GEO accessions
        """

        valid_accessions = list(PreprocessingRunner._get_series_probesets(accessions).keys())

        return valid_accessions

    @staticmethod
    def _get_valid_gdc_projects(project_names):
        """Returns a list of valid GDC projects from project_names"""

        with open("json/gdc_projects.json", encoding="utf-8") as gdc_projects_json:
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

    @staticmethod
    def _get_series_probesets(accessions):
        """Returns a dict of accessions with their respective probesets"""

        with open("json/series_platforms.json", encoding="utf-8") as series_platforms_json:
            series_platforms = json.load(series_platforms_json)

        with open("json/platform_probesets.json", encoding="utf-8") as platform_probesets_json:
            platform_probesets = json.load(platform_probesets_json)

        supported_probesets = set()
        for probemap_fname in [Path(f) for f in os.listdir("probe_maps")]:
            fname_base = probemap_fname.stem
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

    @staticmethod
    def _map_probes(counts_path, species):
        """Maps probes to gene symbols"""

        counts_fname = counts_path.split("/")[-1]
        accession_id = counts_fname.split("_")[0].lower()
        counts_df = pd.read_csv(counts_path, sep="\t")

        probeset = PreprocessingRunner._get_series_probesets([accession_id])[accession_id]
        if not probeset:
            return ""
        counts_df.rename(columns={"probe": probeset}, inplace=True)
        probeset_map_path = f"probe_maps/{species}_{probeset}.tsv"
        if not os.path.exists(probeset_map_path):
            print(f"No probe map found: {probeset_map_path}")
        else:
            print(f"Probe map found: {probeset_map_path}")
        probe_map = pd.read_csv(probeset_map_path, sep="\t",
                                dtype={"hgnc_symbol": object, probeset: object})

        counts_df[probeset] = counts_df[probeset].apply(lambda x: str(x).rstrip("_at"))
        probe_map[probeset] = probe_map[probeset].apply(lambda x: str(x).rstrip("_at"))

        counts_df = pd.merge(counts_df, probe_map, on=probeset, how="outer")

        cols = list(counts_df.columns)
        cols = [cols[-1]] + cols[1:-1]
        counts_df = counts_df[cols]
        counts_df.dropna(inplace=True)

        return counts_df


    def prep_geo_counts(self):
        """Prepares unmapped expression matrices from GEO"""

        counts_paths = self._get_unmapped_counts_paths()
        counts_unmapped = set()
        for counts_path in counts_paths:
            counts_df = PreprocessingRunner._map_probes(counts_path, "hsapiens")
            if not isinstance(counts_df, pd.DataFrame) or len(counts_df) == 0:
                counts_unmapped.add(counts_path)
                continue
            new_counts_path = counts_path.replace("_unmapped.tsv", "_processed.tsv")
            counts_df.to_csv(new_counts_path, sep="\t", index=False)

        # Cleanup unmapped counts
        for counts_path in counts_paths:
            os.remove(counts_path)

        return counts_unmapped

    @staticmethod
    def _zip_preprocessed_data(data_dir):
        """Zips counts and coldata files and returns filename"""

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
                if "counts" in file or "coldata" in file:
                    file_category = "counts" if "counts" in file else "coldata"
                    matched_category = "coldata" if "counts" in file else "counts"
                    matched_file = file.replace(file_category, matched_category)
                    if os.path.exists(file) and os.path.exists(matched_file):
                        print(f"zipping file: {file}")
                        zip_file.write(file)

        # clean up files
        for fpath in files_to_zip:
            if os.path.exists(fpath):
                os.remove(fpath)

        os.chdir(old_wd)

        return zip_fname
