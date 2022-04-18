# Name: Ai Vu

import unittest
import yaml
from webapp import newHelpers
import os


def read_yaml(filename):
    config_file = open(filename, 'r')
    contents = yaml.safe_load(config_file)
    config_file.close()
    return newHelpers.validate_parameters(contents)


class TestParam(unittest.TestCase):
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
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Missing parameter(s): min_expr \n'
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
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Missing parameter(s): min_prop \n'
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
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Missing parameter(s): adj_method \n'
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
            expected = 'Missing parameter(s): condition \n'
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
            expected = 'Missing parameter(s): padj_thresh \n'
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
            expected = 'use_qual_weights must be a bool.\n'
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
            expected = 'Missing parameter(s): contrast_level \n'
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
            expected = 'Missing parameter(s): reference_level \n'
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
            expected = '"padj_thresh" must be a number in range [0,1]\n'
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
            expected = '"min_prop" must be a number in range [0,1]\n'
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
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = '"contrast_level" supported by Microarray analysis is either "normal" or "endometriosis". ' \
                       'You entered "disease"\n'
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
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = 'Currently, program only support condition named "condition". ' \
                       'You entered "I\'m healthy enough"\n'
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
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            # print('result', read_yaml(self.file))
            expected = '"reference_level" supported by Microarray analysis is either "normal" or "endometriosis".' \
                       ' You entered "cancer"\n'
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
        with open(self.file, 'w') as outfile:
            yaml.dump(param, outfile, default_flow_style=False)
            expected = '"contrast_level" supported by RNA-sequence analysis is ' \
                       'either "tumor" or "healthy". You entered "endometriosis"\n' \
                       '"reference_level" supported by RNA-sequence analysis is ' \
                       'either "tumor" or "healthy". You entered "cancer"\n'
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
            expected = '"reference_level" and "contrast_level" cannot refer to the same value.\n'
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
