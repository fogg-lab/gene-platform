# Name: Ai Vu
# Description: Write test for validating input files (tsv) for microarray and rna-seq count and col data

import unittest
from webapp import validate_input_files as valid
from csv import writer
import os


def append_list_as_row(file_name, header, data):
    with open(file_name, 'w', encoding='UTF8', newline='') as write_obj:
        csv_writer = writer(write_obj, delimiter="\t")
        csv_writer.writerow(header)
        for row in data:
            csv_writer.writerow(row)


class TestInput(unittest.TestCase):
    """test check_coldata_matches_counts function which is located in the helpers.py"""

    def setUp(self):
        self.col = 'rna_seq_ucec_coldata.tsv'
        self.count = 'rna_seq_ucec_counts.tsv'

    def tearDown(self):
        try:
            os.remove(self.col)
            os.remove(self.count)
        except FileNotFoundError as fnf:
            print(fnf)

    def test1(self):
        self.col = "rna_col.tsv"
        self.count = "rna_count.tsv"
        header_col = ['sample_name', 'condition', 'data_source']
        field_list_col = [
                          ['GTEX-T6MO-1526-SM-4DM57', 'healthy', 'GTEx'],
                          ['TCGA-FI-A3PX-01A-11R-A22K-07', 'tumor', 'TCGA']]
        header_count = ['Hugo_Symbol', 'Entrez_Gene_Id', 'GTEX-T6MO-1526-SM-4DM57', 'TCGA-FI-A3PX-01A-11R-A22K-07']
        field_list_count = [
                            ['FAM208A', 23272, 2523.95, 183.0], ['RADIL', 55698, 320.0, 10.98273265687188]]
        validate = valid.FileValidation()
        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        self.assertEqual(validate.validate_file(self.col, self.count)[1], "")

    def test2(self):
        self.col = 'rna_col.tsv'
        self.count = 'rna_count.tsv'
        col_head = ['sample_name', 'condition', 'data_source']
        field_list_col = [
                          ['GTEX-T6MO-1526-SM-4DM57', 'healthy', 'GTEx'],
                          ['TCGA-FI-A3PX-01A-11R-A22K-07', 'tumor', 'TCGA']]

        count_head = ['Hugo_Symbol', 'Entrez_Gene_Id', 'GTEX-T6MO-1526-SM-4DM57', 'GTEX-11P81-1626-SM-5BC52']
        field_list_count = [['FAM208A', 23272, 2523.95, 183.0], ['RADIL', 55698, 320.0, 10.98273265687188]]
        append_list_as_row(self.col, col_head, field_list_col)
        append_list_as_row(self.count, count_head, field_list_count)

        validate = valid.FileValidation()

        expected = "GTEX-11P81-1626-SM-5BC52 presents in rna-seq col data file but not in rna-seq count data file\n"
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expected)

    def test3(self):
        self.col = 'rna_col.tsv'
        self.count = 'rna_count.tsv'
        header_col = ['sample_name', 'condition', 'data_source']
        field_list_col = [
                          ['GTEX-T6MO-1526-SM-4DM57', 'healthy', 'GTEx'],
                          [1, 'tumor', 'TCGA']]
        header_count = ['Hugo_Symbol', 'Entrez_Gene_Id', 'GTEX-T6MO-1526-SM-4DM57', 1]
        field_list_count = [
            ['FAM208A', 23272, 2523.95, 183.0], ['RADIL', 55698, 320.0, 10.98273265687188]]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        self.assertEqual(validate.validate_file(self.col, self.count)[1], "")

    def test4(self):
        self.col = 'rna_col.tsv'
        self.count = 'rna_count.tsv'
        header_col = ['sample_name', 'condition', 'data_source']
        field_list_col = [
            ['GTEX-T6MO-1526-SM-4DM57', 'healthy', 'GTEx'],
            [1, 'tumor', 'TCGA']]
        header_count = ['Hugo_Symbol', 'Entrez_Gene_Id', 'GTEX-T6MO-1526-SM-4DM57', 1]
        field_list_count = [
            ['FAM208A', 23272.12, 2523.95, 183.0], ['RADIL', 55698, 320.0, 10.98273265687188]]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = "The column Entrez_Gene_Id in rna_seq_count file has data type of float64, " \
                 "which should have been int64\n"
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)

    def test5(self):
        self.col = 'rna_col.tsv'
        self.count = 'rna_count.tsv'
        header_col = ['sample_name', 'condition', 'data_source']
        field_list_col = [
            ['GTEX-T6MO-1526-SM-4DM57', 'healthy', 'GTEx'],
            [1, 'tumor', 'TCGA']]
        header_count = ['Hugo_Symbol', 'Entrez_Gene_Id', 'GTEX-T6MO-1526-SM-4DM57', 1]
        field_list_count = [
            ['FAM208A', 23272, 'hi', 183.0], ['RADIL', 55698, 320.0, 10.98273265687188]]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = "The column GTEX-T6MO-1526-SM-4DM57 in rna_seq_count file has " \
                 "data type of object, which should have been float64\n"
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)


    def test6(self):
        self.col = 'rna_col.tsv'
        self.count = 'rna_count.tsv'
        header_col = ['sample_name', 'condition', 'data_source']
        field_list_col = [
            ['GTEX-T6MO-1526-SM-4DM57', 'healthy', 'GTEx'],
            [1, 'tumor', 'TCGA']]
        header_count = ['Hugo_Symbol', 'Entrez_Gene_Id', 'GTEX-T6MO-1526-SM-4DM57', 1]
        field_list_count = [
            ['FAM208A', 23272, 249, 183.0], ['RADIL', 55698, 320.0, 10.98273265687188]]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = ""
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)

    def test7(self):
        self.col = 'rna_col.tsv'
        self.count = 'rna_count.tsv'
        header_col = ['sample_name', 'condition', 'data_source']
        field_list_col = [
            ['GTEX-T6MO-1526-SM-4DM57', 'healthy', 'GTEx'],
            [1, 'tumor', 'TCGA']]
        header_count = ['Hugo_Symbol', 'Entrez_Gene_Id', 'GTEX-T6MO-1526-SM-4DM57', 1]
        field_list_count = [
            ['FAM208A', 23272, 183.0], ['RADIL', 55698, 320.0, 10.98273265687188]]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = "1 missing value(s) in your count data file\n"
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)

    def test8(self):
        self.col = 'micro_col.tsv'
        self.count = 'micro_count.tsv'
        header_col = ['sample_name', 'series', 'title', 'condition', 'endometriosis_stage', 'phase', 'tissue', 'batch']
        field_list_col = [
            ['GSM109815', 'GSE4888', 'endometrium M165', 'normal', 'none', 'proliferative', 'uterus', 1],
            ['GSM150190', 'GSE6364', 'PE_D_26A', 'endometriosis', 'moderate_severe', 'proliferative', 'uterus',1]]
        header_count = ['symbol', 'GSM109815', 'GSM150190']
        field_list_count = [
            ['A1BG-AS1', 6.354643206196131, 7.02545678991068], ['A2ML1', 10.70514655195897, 11.434625316952236]]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = ""
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)

    def test9(self):
        self.col = 'micro_col.tsv'
        self.count = 'micro_count.tsv'
        header_col = ['sample_name', 'series', 'title', 'condition', 'endometriosis_stage', 'phase', 'tissue',
                      'batch']
        field_list_col = [
            ['GSM109815', 'GSE4888', 'endometrium M165', 'normal', 'none', 'proliferative', 'uterus', 1],
            ['GSM150190', 'GSE6364', 'PE_D_26A', 'endometriosis', 'moderate_severe', 'proliferative', 'uterus', 1]]
        header_count = ['symbol', 'GSM109815', 'GSM150190']
        field_list_count = [
            ['A1BG-AS1', 7.02545678991068], ['A2ML1', 10.70514655195897, 11.434625316952236]]
        validate = valid.FileValidation()
        # validate = FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = "1 missing value(s) in your count data file\n"
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)

    def test10(self):
        self.col = 'micro_col.tsv'
        self.count = 'micro_count.tsv'
        header_col = ['sample_name', 'series', 'title', 'condition', 'endometriosis_stage', 'phase', 'tissue',
                      'batch']
        field_list_col = [
            ['GSM109815', 'GSE4888', 'endometrium M165', 'normal', 'none', 'proliferative', 'uterus', 1],
            ['GSM150190', 'GSE6364', 'PE_D_26A', 'endometriosis', 'moderate_severe', 'proliferative', 'uterus', 1]]
        header_count = ['symbol', 'GSM109815', 'GSM150190']
        field_list_count = [
            [1, 6, 7.02545678991068], ['A2ML1', 10.70514655195897, 11.434625316952236]]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = ""
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)


    def test11(self):
        self.col = 'micro_col.tsv'
        self.count = 'micro_count.tsv'
        header_col = ['sample_name', 'series', 'title', 'endometriosis_stage', 'phase', 'tissue',
                      'batch']
        field_list_col = [
            ['GSM109815', 'GSE4888', 'endometrium M165', 'normal', 'none', 'proliferative', 'uterus', 1],
            ['GSM150190', 'GSE6364', 'PE_D_26A', 'endometriosis', 'moderate_severe', 'proliferative', 'uterus', 1]]
        header_count = ['symbol', 'GSM109815', 'GSM150190']
        field_list_count = [
            [1, 6, 7.02545678991068], ['A2ML1', 10.70514655195897, 11.434625316952236]]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = "column names of your submitted cold data file:\n" \
                 "['sample_name', 'series', 'title', 'endometriosis_stage', 'phase', 'tissue', 'batch']\n" \
                 "column names of expected cold file for microarray:\n" \
                 "['sample_name', 'series', 'title', 'condition', 'endometriosis_stage', 'phase', 'tissue', 'batch']\n" \
                 "column names of expected cold file for rna_sequence:\n" \
                 "['sample_name', 'condition', 'data_source']\n"
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)

    def test12(self):
        self.col = 'micro_col.tsv'
        self.count = 'micro_count.tsv'
        header_col = ['sample_name', 'series', 'title', 'condition', 'endometriosis_stage', 'phase', 'tissue',
                      'batch']
        field_list_col = [
            ['GSM109815', 'GSE4888', 'endometrium M165', 'normal', 'none', 'proliferative', 'uterus', 1],
            ['GSM150190', 'GSE6364', 'PE_D_26A', 'endometriosis', 'moderate_severe', 'proliferative', 'uterus', 1]]
        header_count = ['symbol', 'GSM109815', 'GSM150190']
        field_list_count = [
            ["", 7.02545678991068], ['A2ML1', 10.70514655195897, 11.434625316952236]]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = "2 missing value(s) in your count data file\n"
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)

    def test13(self):
        self.col = 'micro_col.tsv'
        self.count = 'micro_count.tsv'
        header_col = ['sample_name', 'series', 'title', 'condition', 'endometriosis_stage', 'phase', 'tissue',
                      'batch']
        field_list_col = [

        ]
        header_count = ['symbol', 'GSM109815', 'GSM150190']
        field_list_count = [
            ["", 7.02545678991068], ['A2ML1', 10.70514655195897, 11.434625316952236]]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = "2 missing value(s) in your count data file\n" \
                 "Microarray col data file has 0 row of data\n"
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)

    def test14(self):
        self.col = 'micro_col.tsv'
        self.count = 'micro_count.tsv'
        header_col = ['sample_name', 'series', 'title', 'condition', 'endometriosis_stage', 'phase', 'tissue',
                      'batch']
        field_list_col = [

        ]
        header_count = ['symbol', 'GSM109815', 'GSM150190']
        field_list_count = [

        ]
        validate = valid.FileValidation()

        append_list_as_row(self.col, header_col, field_list_col)
        append_list_as_row(self.count, header_count, field_list_count)
        expect = "Microarray col data file has 0 row of data\n" \
                 "Microarray count data file has 0 row of data\n"
        self.assertEqual(validate.validate_file(self.col, self.count)[1], expect)


if __name__ == "__main__":
    unittest.main()