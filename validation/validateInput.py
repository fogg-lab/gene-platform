# Name: Ai Vu
# Description: This program is used to validata the 5 input files: microarray coldata, microarray countdata,
# rna-seq coldata, rna-seq countdata, and filter gene list. To validate, users create a FileValidation object
# and use the function validate_file with the filename as an argument. There are 3 possible result: True if
# file is valid, False if file is invalid, and Pending if there is not enough information to conclude. Pending
# does not mean file is invalid, but it might be due to the FileValidation needs to wait for the other file (either
# coldata or count data of the same analysis) to complete its task.

import pandas as pd
import numpy as np
import os
from cerberus import Validator
import yaml

UPLOAD_FOLDER = r"C:\Users\trungn\PycharmProjects\DGEAP\validation"


class FileSystem:
    """holds valid formats of col_data, counts, and filter gene file"""

    def __init__(self):
        self.micro_array_col_data_dict = {'sample_name': "string", 'series': "string", 'title': "string",
                                          'condition': "string", 'endometriosis_stage': "string",
                                          'phase': "string", 'tissue': "string", 'batch': int}
        self.rna_seq_col_dict = {'sample_name': "string", 'condition': "string", 'data_source': "string"}
        self.filter_gene_dict = {0: "string"}
        self.micro_array_count_data_dict = {'symbol': "string"}
        self.rna_seq_count_data_dict = {'Hugo_Symbol': 'string', 'Entrez_Gene_Id': int}

        self._filter_gene_file = {'type': '.txt', 'format': self.filter_gene_dict}

        self._micro_array_col_file = {'type': '.tsv', 'format': self.micro_array_col_data_dict}

        self._micro_array_count_file = {'type': '.tsv', 'format': self.micro_array_count_data_dict}

        self._rna_seq_col_file = {'type': '.tsv', 'format': self.rna_seq_col_dict}

        self._rna_seq_count_file = {'type': '.tsv', 'format': self.rna_seq_count_data_dict}
        self._microarray_scheme = {
            'min_expr': {
                'type': 'integer',
                'required': True
            },
            'min_prop': {
                'type': 'float', 'min': 0.0, 'max': 1.0,
                'required': True
            },
            'adj_method': {
                'type': 'string', 'allowed': ['BH'],
                'required': True
            },
            'condition': {
                'type': 'string', 'allowed': ['condition'],
                'required': True
            },
            'use_weight': {
                'type': 'boolean',
                'required': True
            },
            'contrast_level': {
                'type': 'string', 'allowed': ['endometriosis', 'normal'],
                'required': True
            },
            'reference_level': {
                'type': 'string', 'allowed': ['endometriosis', 'normal'],
                'required': True
            },
        }

        self._rna_seq_schema = {
            'min_expr': {
                'type': 'integer',
                'required': True
            },
            'min_prop': {
                'type': 'float', 'min': 0.0, 'max': 1.0,
                'required': True
            },
            'padj_thresh': {
                'type': 'float', 'min': 0.0, 'max': 1.0,
                'required': True
            },
            'adj_method': {
                'type': 'string', 'allowed': ['BH'],
                'required': True
            },
            'condition': {
                'type': 'string', 'allowed': ['condition'],
                'required': True
            },
            'contrast_level': {
                'type': 'string', 'allowed': ['tumor', 'healthy'],
                'required': True
            },
            'reference_level': {
                'type': 'string', 'allowed': ['tumor', 'healthy'],
                'required': True
            }
        }

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

    def get_file_format_filter_gene(self):
        """get the valid format of the filter gene file"""
        return self._filter_gene_file

    def get_config_schema_microarray(self):
        """get the valid format of the config file for microarray's parameters"""
        return self._microarray_scheme

    def get_config_schema_rna_seq(self):
        """get the valid format of the config file for microarray's parameters"""
        return self._rna_seq_schema


class FileValidation:
    """has the FileSystem and methods to validate input files: col-data, count data, and filter gene files"""

    def __init__(self):

        self._validate_status = False  # assume that the file has invalid format
        self._has_validated_micro_array_col_file = False  # holds the validated micro_array filename
        self._has_validated_rna_seq_col_file = False  # # holds the validated rna_seq filename
        self.file_obj = FileSystem()
        self._filename = ['filter gene', 'rna seguence coldata', 'rna sequence_count',
                          'microarray coldata', 'microarray_count', 'config_microarray']

    def set_validation(self, state):
        """set the new validation status for a file"""
        self._validate_status = state

    def get_validation(self):
        """get the validation status of a file"""
        return self._validate_status

    def error_message(self, df, filename):
        """print valid format and format of the file when there is type error"""
        print("the file has type error. Below is the data types in your files")
        print(df.dtypes)
        print("The valid data types format is shown as below:")
        if filename == self._filename[0]:
            print(pd.DataFrame.from_dict(self.file_obj.get_file_format_filter_gene()))
        if filename == self._filename[1]:
            print(pd.DataFrame.from_dict(self.file_obj.get_file_format_rna_seq_col()))
        if filename == self._filename[2]:
            print(pd.DataFrame.from_dict(self.file_obj.get_file_format_rna_seq_count()))
        if filename == self._filename[3]:
            print(pd.DataFrame.from_dict(self.file_obj.get_file_format_micro_array_col()))
        if filename == self._filename[4]:
            print(pd.DataFrame.from_dict(self.file_obj.get_file_format_micro_array_count()))

    @staticmethod
    def pending_message(filename):
        return "pending: waiting for the " + filename + " to finish validation"

    @staticmethod
    def convert_string_to_dict(column_list, data_type):
        """convert a list into dictionary with values be data type"""
        string_data_type_dict = {}
        for gene in column_list:
            string_data_type_dict[gene] = data_type
        return string_data_type_dict

    def validate_config(self, config_name):
        file_path = os.path.join(UPLOAD_FOLDER, config_name)
        try:
            with open(file_path, 'r') as stream:
                try:
                    config = yaml.safe_load(stream)
                    if config_name.endswith('microarray.yml'):
                        v = Validator(self.file_obj.get_config_schema_microarray())
                        result = v.validate(config, self.file_obj.get_config_schema_microarray())
                        if result is True:
                            self.set_validation(True)
                        else:
                            self.set_validation(False)
                            print(v.errors)

                    if config_name.endswith('rna_seq.yml'):
                        v = Validator(self.file_obj.get_config_schema_rna_seq())
                        result = v.validate(config, self.file_obj.get_config_schema_rna_seq())
                        if result is True:
                            self.set_validation(True)
                        else:
                            self.set_validation(False)
                            print(v.errors)
                except yaml.YAMLError as exception:
                    raise exception

        except IOError:
            fatal_error("Could not open the configuration file name {}".format(config_name))

        finally:
            return self.get_validation()

    def validate_file(self, filename):
        """validate the following files: col, count, and filter gene. Returns true
         if file is valid; otherwise, returns false"""

        file_path = os.path.join(UPLOAD_FOLDER, filename)
        df = pd.DataFrame()
        name = ''
        try:
            if filename.endswith(".txt"):
                name = self._filename[0]
                df = pd.read_csv(file_path, sep=" ", header=None,
                                 dtype=self.file_obj.get_file_format_filter_gene()['format'])
                self.set_validation(True)

            elif filename.endswith(".tsv"):
                df = pd.read_csv(file_path, sep='\t')
                if len(df.columns) == len(self.file_obj.get_file_format_micro_array_col()['format']):
                    name = self._filename[3]
                    df = pd.read_csv(file_path, sep='\t',
                                     dtype=self.file_obj.get_file_format_micro_array_col()['format'])
                    self.file_obj.merge_dicts(
                        self.convert_string_to_dict(df['sample_name'].tolist(), np.float),
                        self.file_obj.micro_array_count_data_dict)
                    self._has_validated_micro_array_col_file = True
                    # print("micro_count", self.file_obj.get_file_format_micro_array_count())
                    self.set_validation(True)

                elif len(df.columns) == len(self.file_obj.get_file_format_rna_seq_col()['format']):
                    name = self._filename[1]
                    df = pd.read_csv(file_path, sep='\t', dtype=self.file_obj.get_file_format_rna_seq_col()['format'])
                    self.file_obj.merge_dicts(
                        self.convert_string_to_dict(df['sample_name'].tolist(), np.float),
                        self.file_obj.rna_seq_count_data_dict)
                    # print("rna_count", self.file_obj.get_file_format_rna_seq_count())
                    self._has_validated_rna_seq_col_file = True
                    self.set_validation(True)

                elif df.columns.values.tolist()[0] in self.file_obj.get_file_format_micro_array_count()['format']:
                    if self._has_validated_micro_array_col_file is True:
                        name = self._filename[4]
                        df = pd.read_csv(file_path, sep='\t',
                                         dtype=self.file_obj.get_file_format_micro_array_count()['format'])
                        self.set_validation(True)
                    else:
                        self.set_validation(self.pending_message(self._filename[3]))
                elif df.columns.values.tolist()[0] in self.file_obj.get_file_format_rna_seq_count()['format']:
                    if self._has_validated_rna_seq_col_file is True:
                        name = self._filename[2]
                        df = pd.read_csv(file_path, sep='\t',
                                         dtype=self.file_obj.get_file_format_rna_seq_count()['format'])
                        self.set_validation(True)
                    else:
                        self.set_validation(self.pending_message(self._filename[1]))

        except TypeError:
            self.error_message(df, name)
            self.set_validation(False)

        finally:
            return self.get_validation()


# def main():
#     # this is to show examples of how to use the program to check for file validation
#     valid = FileValidation()
#     config1 = 'config_microarray.yml'
#     print(valid.validate_config(config1))
#
#     config2 = 'config_rna_seq.yml'
#     print(valid.validate_config(config2))
#
#     file_name_1 = "filter_gene.txt"
#     print(file_name_1)
#     print(valid.validate_file(file_name_1))
#
#     print("2= expect pend=========================")
#     file_name_4 = 'microarray_endometriosis_counts.tsv'
#     print(file_name_4)
#     print(valid.validate_file(file_name_4))
#
#     print("3= expect pend=========================")
#     file_name_5 = 'rna_seq_ucec_counts.tsv'
#     print(file_name_5)
#     print(valid.validate_file(file_name_5))
#
#     print("4=true=========================")
#     file_name_2 = 'microarray_endometriosis_coldata.tsv'
#     print(file_name_2)
#     print(valid.validate_file(file_name_2))
#
#     print("6=true=========================")
#     file_name_3 = 'rna_seq_ucec_coldata.tsv'
#     print(file_name_3)
#     print(valid.validate_file(file_name_3))
#
#     # print("5=true=========================")
#     # file_name_4 = 'microarray_endometriosis_counts.tsv'
#     # print(file_name_4)
#     # print(valid.validate_file(file_name_4))
#     # print("==========================")
#
#     # print("7=true=========================")
#     # file_name_5 = 'rna_seq_ucec_counts.tsv'
#     # print(file_name_5)
#     # print(valid.validate_file(file_name_5))
#
#
# if __name__ == '__main__':
#     main()
