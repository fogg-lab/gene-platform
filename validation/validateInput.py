import pandas as pd
import numpy as np
import os

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


class FileValidation:
    """has the FileSystem and methods to validate input files: col-data, count data, and filter gene files"""
    def __init__(self):

        self._validate_status = False  # assume that the file has invalid format
        self._has_validated_micro_array_col_file = False  # holds the validated micro_array filename
        self._has_validated_rna_seq_col_file = False  # # holds the validated rna_seq filename
        self.file_obj = FileSystem()
        self._filename = ['filter gene', 'rna seguence coldata', 'rna sequence_count',
                          'microarray coldata', 'microarray_count']

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
        dict_1 = {}
        for gene in column_list:
            dict_1[gene] = data_type
        return dict_1

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


def main():
    valid = FileValidation()
    file_name_1 = "filter_gene.txt"
    print(file_name_1)
    print(valid.validate_file(file_name_1))

    print("2= expect pend=========================")
    file_name_4 = 'microarray_endometriosis_counts.tsv'
    print(file_name_4)
    print(valid.validate_file(file_name_4))

    print("3= expect pend=========================")
    file_name_5 = 'rna_seq_ucec_counts.tsv'
    print(file_name_5)
    print(valid.validate_file(file_name_5))

    print("4=true=========================")
    file_name_2 = 'microarray_endometriosis_coldata.tsv'
    print(file_name_2)
    print(valid.validate_file(file_name_2))

    print("5=true=========================")
    file_name_4 = 'microarray_endometriosis_counts.tsv'
    print(file_name_4)
    print(valid.validate_file(file_name_4))
    print("==========================")

    print("6=true=========================")
    file_name_3 = 'rna_seq_ucec_coldata.tsv'
    print(file_name_3)
    print(valid.validate_file(file_name_3))

    print("7=true=========================")
    file_name_5 = 'rna_seq_ucec_counts.tsv'
    print(file_name_5)
    print(valid.validate_file(file_name_5))


if __name__ == '__main__':
    main()

