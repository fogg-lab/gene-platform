# Name: Ai Vu
# Description: Write test to validate parameters (yaml file) for microarray and rna-seq
# and the check_coldata_matches_counts helper function.

import unittest
import yaml
from webapp import helpers
import csv


def read_yaml(filename):
    config_file = open(filename, 'r')
    contents = yaml.safe_load(config_file)
    config_file.close()
    return helpers.validate_parameters(contents)


def get_inputs_files(filename):
    print(filename)
    with open(filename, newline='') as inp:
        input = csv.reader(inp, delimiter="\t")
        input_as_list = list(input)
        return input_as_list


def test_input_files(col_file, count_file):
    """tests all input .tsv files"""
    col_data_list = get_inputs_files(col_file)
    counts_list = get_inputs_files(count_file)

    confirmation_message = helpers.check_coldata_matches_counts(
        copy.deepcopy(counts_list[0]), copy.deepcopy(col_data_list))
    return confirmation_message


class TestParam(unittest.TestCase):
    """test the check_factor_level functions and validate_parameters function
    which are located in helpers.py"""
    def tearDown(self):
        try:
            os.remove(self.file)
        except FileNotFoundError as fnf:
            print(fnf)

    def test1(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.25,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.05,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="normal"
        )
        self.assertEqual(helpers.check_factor_levels(param, get_inputs_files("microarray_endometriosis_coldata.tsv")), "")
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = ''
            self.assertEqual(read_yaml(self.file), expected)

    def test2(self):
        self.file = 'config.yaml'
        param = dict(
            min_prop=0.25,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.05,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="normal"
        )
        self.assertEqual(helpers.check_factor_levels(param, get_inputs_files("microarray_endometriosis_coldata.tsv")), "")
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Missing parameter: min_expr\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test3(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            # min_prop=0.25,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.05,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="normal"
        )
        self.assertEqual(helpers.check_factor_levels(param, get_inputs_files("microarray_endometriosis_coldata.tsv")), "")
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Missing parameter: min_prop\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test4(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.25,
            # adj_method="BH",
            condition="condition",
            padj_thresh=0.05,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="normal"
        )
        self.assertEqual(helpers.check_factor_levels(param, get_inputs_files("microarray_endometriosis_coldata.tsv")), "")
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Missing parameter: adj_method\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test5(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.25,
            adj_method="BH",
            # condition="condition",
            padj_thresh=0.05,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="normal"
        )
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Missing parameter: condition\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test6(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.25,
            adj_method="BH",
            condition="condition",
            # padj_thresh=0.05,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="normal"
        )
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Missing parameter: padj_thresh\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test7(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.25,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.05,
            use_qual_weights=0.5,
            contrast_level="endometriosis",
            reference_level="normal"
        )
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'use_qual_weights must be either True or False.\n'
            self.assertEqual(read_yaml(self.file), expected)  # empty. logic?

    def test8(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.25,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.05,
            use_qual_weights=True,
            # contrast_level="endometriosis",
            reference_level="normal"
        )
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Missing parameter: contrast_level\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test9(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.25,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.05,
            use_qual_weights=True,
            contrast_level="endometriosis",
            # reference_level="normal"
        )
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Missing parameter: reference_level\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test10(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.25,
            adj_method="BH",
            condition="condition",
            padj_thresh=1.05,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="normal"
        )
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            print('result', read_yaml(self.file))
            expected = '"padj_thresh" must be between 0 and 1\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test11(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=1.01,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.5,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="normal"
        )
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = '"min_prop" must be between 0 and 1\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test12(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.78,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.5,
            use_qual_weights=True,
            # contrast_level="endometriosis",
            contrast_level='disease',
            reference_level="normal"
        )
        print('result', helpers.check_factor_levels(param, get_inputs_files("microarray_endometriosis_coldata.tsv")))
        self.assertEqual(helpers.check_factor_levels(param, get_inputs_files("microarray_endometriosis_coldata.tsv")),
                         "Unknown contrast level 'disease'\n")

        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = ''
            self.assertEqual(read_yaml(self.file), expected)

    def test13(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.78,
            adj_method="BH",
            condition="I'm healthy enough",
            padj_thresh=0.5,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="normal"
        )
        self.assertEqual(helpers.check_factor_levels(param, get_inputs_files("microarray_endometriosis_coldata.tsv")),
                         "Condition 'I'm healthy enough' not present in coldata's column name\n")
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = ''
            self.assertEqual(read_yaml(self.file), expected)

    def test14(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.78,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.5,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="cancer"
        )
        self.assertEqual(helpers.check_factor_levels(param, get_inputs_files("microarray_endometriosis_coldata.tsv")),
                         "Unknown reference level 'cancer'\n")

        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = ''
            self.assertEqual(read_yaml(self.file), expected)

    def test15(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.78,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.5,
            # use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="cancer"
        )
        self.assertEqual(helpers.check_factor_levels(param, get_inputs_files("microarray_endometriosis_coldata.tsv")),
                         "Unknown reference level 'cancer'\n")

        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = ''
            self.assertEqual(read_yaml(self.file), expected)

    def test16(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.78,
            adj_method="BH",
            condition="condition",
            padj_thresh=0.5,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="endometriosis"
        )
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'reference_level and contrast_level cannot be the same.\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test18(self):
        self.file = 'config.yaml'
        param = dict(
            min_expr=5.6438,
            min_prop=0.78,
            adj_method="I don't know",
            condition="condition",
            padj_thresh=0.5,
            use_qual_weights=True,
            contrast_level="endometriosis",
            reference_level="endometriosis"
        )
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            print('result', read_yaml(self.file))
            expected = 'Currently, program only support adj_method named "BH". You entered "I don\'t know".\n' \
                       '"reference_level" and "contrast_level" cannot refer to the same value.\n'
            self.assertEqual(read_yaml(self.file), expected)

    def test18(self):
        self.file = 'config.yaml'
        try:
            param = dict(
                min_expr,
                min_prop=0.78,
                adj_method="BH",
                condition="condition",
                padj_thresh=0.5,
                use_qual_weights=True,
                contrast_level="endometriosis",
                reference_level="normal"
            )
            with open(self.file, 'w') as outfile:
                yaml.dump(param, outfile, default_flow_style=False)
                expected = 'Missing values for parameter(s): min_expr\n'
                self.assertEqual(read_yaml(self.file), expected)

        except NameError as e:
            expected = "name 'min_expr' is not defined"
            self.assertEqual(str(e), expected)

    def test19(self):
        self.file = 'config.yaml'
        try:
            param = dict(
                0.7,
                min_prop=0.78,
                adj_method="BH",
                condition="condition",
                padj_thresh=0.5,
                use_qual_weights=True,
                contrast_level="endometriosis",
                reference_level="normal"
            )
            with open(self.file, 'w') as outfile:
                yaml.dump(param, outfile, default_flow_style=False)
                expected = 'Missing values for parameter(s): min_expr\n'
                self.assertEqual(read_yaml(self.file), expected)
        except TypeError as e:
            expected = "'float' object is not iterable"
            self.assertEqual(str(e), expected)


if __name__ == "__main__":
    unittest.main()
