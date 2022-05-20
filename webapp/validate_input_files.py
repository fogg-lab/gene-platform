# Name: Ai Vu

import pandas as pd
import numpy as np
import os

# UPLOAD_FOLDER = r"C:\Users\trungn\PycharmProjects\DGEAP"


class FileSystem:
    """holds valid formats of col_data, counts, and filter gene file"""

    def __init__(self):
        self.micro_array_col_data_dict = {'sample_name': "string", 'series': "string", 'title': "string",
                                          'condition': "string", 'endometriosis_stage': "string",
                                          'phase': "string", 'tissue': "string", 'batch': int}
        self.rna_seq_col_dict = {'sample_name': "string", 'condition': "string", 'data_source': "string"}
        self.micro_array_count_data_dict = {'symbol': "string"}
        self.rna_seq_count_data_dict = {'Hugo_Symbol': 'object', 'Entrez_Gene_Id': int}
        self._micro_array_col_file = {'type': '.tsv', 'column_count': len(self.micro_array_col_data_dict),
                                      'format': self.micro_array_col_data_dict}

        self._micro_array_count_file = {'type': '.tsv', 'column_count': len(self.micro_array_count_data_dict),
                                        'format': self.micro_array_count_data_dict}

        self._rna_seq_col_file = {'type': '.tsv', 'column_count': len(self.rna_seq_col_dict),
                                  'format': self.rna_seq_col_dict}

        self._rna_seq_count_file = {'type': '.tsv', 'column_count': len(self.rna_seq_count_data_dict),
                                    'format': self.rna_seq_count_data_dict}

    @staticmethod
    def merge_dicts(dict_1, final_dict):
        """combine two dicts into one by appending"""
        final_dict.update(dict_1)

    def get_file_format_micro_array_col(self):
        """get the valid format for micro_array col_data file"""
        return self._micro_array_col_file

    def get_file_format_micro_array_count(self):
        """get the valid format for micro_array count file"""
        return self._micro_array_count_file

    def get_file_format_rna_seq_col(self):
        """get the valid format of rna_seq col_data file"""
        return self._rna_seq_col_file

    def get_file_format_rna_seq_count(self):
        """get the valid format of the rna_seq count file"""
        return self._rna_seq_count_file


class FileValidation:
    """has the FileSystem and methods to validate input files: col-data, count data, and filter gene files"""

    def __init__(self):

        self._validate_status = False  # assume that the file has invalid format
        self._has_validated_micro_array_col_file = False  # holds the validated micro_array filename
        self._has_validated_rna_seq_col_file = False  # # holds the validated rna_seq filename
        self.file_obj = FileSystem()

    def set_validation(self, state):
        """set the new validation status for a file"""
        self._validate_status = state

    def get_validation(self):
        """get the validation status of a file"""
        return self._validate_status

    def error_message(self, df, file):
        """print valid format and format of the file when there is type error"""
        y = {k: str(v[0]) for k, v in pd.DataFrame(df.dtypes).T.to_dict('list').items()}
        if file == 'rna_seq_col':
            x = self.file_obj.get_file_format_rna_seq_col()['format']
        elif file == 'rna_seq_count':
            x = self.file_obj.get_file_format_rna_seq_count()['format']
        elif file == 'microarray_col':
            x = self.file_obj.get_file_format_micro_array_col()['format']
        elif file == 'microarray_count':
            x = self.file_obj.get_file_format_micro_array_count()['format']
        for k in x.keys():
            if x[k] == type(0.0):
                x[k] = 'float64'
            elif x[k] == type(0):
                x[k] = 'int64'
        diff = {k: x[k] for k in x if k in y and x[k] != y[k]}
        msg = ""
        for key, value in diff.items():
            msg += 'The column {} in {} file has data type of {}, which should have been {}\n'.format(key, file, y[key],
                                                                                                      value)
        return msg

    @staticmethod
    def convert_string_to_dict(column_list, expected_datatype):
        """convert a list into dictionary with values be data type"""
        dict_1 = {}
        for gene in column_list:
            dict_1[gene] = expected_datatype
        return dict_1

    def validate_file(self, col_path, count_path):
        """validate the following files: col, count, and filter gene.
        Take file paths of col and count files, returns a tuple (true/false, error message"""

        data_check_col_df = pd.DataFrame()
        data_check_count_df = pd.DataFrame()
        name = ''
        msg = ''
        if not col_path:
            msg = 'there is no coldata file path\n'
        if not count_path:
            msg += 'there is no count data file path\n'
        try:
            col_df = pd.read_csv(col_path, sep='\t')
            count_df = pd.read_csv(count_path, sep='\t')
            missing_count = count_df.isnull().sum().sum()
            missing_col = col_df.isnull().sum().sum()

            if missing_count != 0:
                msg += f'{missing_count} missing value(s) in your count data file\n'
                self.set_validation(False)
            if missing_col !=0:
                msg += f'{missing_col} missing value(s) in your col data file\n'
                self.set_validation(False)

            if col_df.columns.tolist() == list(self.file_obj.get_file_format_micro_array_col()['format'].keys()):
                keys = list(self.file_obj.get_file_format_micro_array_col()['format'].keys())
                name = 'microarray_col'
                data_check_col_df = pd.read_csv(col_path,
                                                dtype=self.file_obj.get_file_format_micro_array_col()['format'],sep='\t')
                if len(data_check_col_df) == 0:
                    msg += f'Microarray col data file has 0 row of data\n'
                    self.set_validation(False)
                self.file_obj.merge_dicts(
                    self.convert_string_to_dict(col_df['sample_name'].tolist(), float),
                    self.file_obj.micro_array_count_data_dict)
                self._has_validated_micro_array_col_file = True
                if count_df.columns.values.tolist()[0] == "symbol":
                    name = 'microarray_count'
                    data_check_count_df = pd.read_csv(count_path,
                                                      dtype=self.file_obj.get_file_format_micro_array_count()['format'], sep='\t')
                    if len(data_check_count_df) == 0:
                        msg += f'Microarray count data file has 0 row of data\n'
                        self.set_validation(False)
                    self.set_validation(True)

            elif col_df.columns.tolist() == list(self.file_obj.get_file_format_rna_seq_col()['format'].keys()):
                name = 'rna_seq_col'

                data_check_col_df = pd.read_csv(col_path, dtype=self.file_obj.get_file_format_rna_seq_col()['format'], sep='\t')
                if len(data_check_col_df) == 0:
                    msg += f'RNA-sequence col data file has 0 row of data\n'
                    return self.get_validation(), msg

                self.file_obj.merge_dicts(
                    self.convert_string_to_dict(col_df['sample_name'].tolist(), float),
                    self.file_obj.rna_seq_count_data_dict)
                self._has_validated_rna_seq_col_file = True
                if count_df.columns.values.tolist()[0] == "Hugo_Symbol":
                    name = 'rna_seq_count'
                    data_check_count_df = pd.read_csv(count_path,
                                                      dtype=self.file_obj.get_file_format_rna_seq_count()['format'], sep="\t")
                    if len(data_check_count_df) == 0:
                        msg += f'RNA-sequence count data file has 0 row of data\n'
                        self.set_validation(False)

                    keys_list = self.file_obj.get_file_format_rna_seq_count()['format'].keys()

                    if keys_list == data_check_count_df.columns.values.tolist():
                        self.set_validation(True)
                    else:
                        for field in data_check_count_df.columns.values.tolist():
                            if field not in keys_list:
                                msg += f"{field} presents in rna-seq count data file but not in " \
                                       f"rna-seq col data file\n"
                                self.set_validation(False)

            else:
                msg += f'column names of your submitted col file:\n{col_df.columns.tolist()}\n'
                msg += f"column names of expected col file for " \
                       f"microarray:\n{list(self.file_obj.get_file_format_micro_array_col()['format'].keys())}\n"
                msg += f"column names of expected col file for " \
                       f"rna_sequence:\n{list(self.file_obj.get_file_format_rna_seq_col()['format'].keys())}\n"

                self.set_validation(False)

        except:
            if len(data_check_col_df) == 0:
                msg += self.error_message(col_df, name)
            if len(data_check_count_df) == 0:
                msg += self.error_message(count_df, name)
            self.set_validation(False)

        finally:
            return self.get_validation(), msg
