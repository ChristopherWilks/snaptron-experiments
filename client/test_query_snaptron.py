#!/usr/bin/env python2.7
import sys
import argparse
import unittest

import query_snaptron

samples_records = None
args = None

def setUpModule():
    global sample_records
    global args
    parser = query_snaptron.create_parser(disable_header=False)
    args = parser.parse_args()
    sample_records = query_snaptron.download_sample_metadata(args)
    pass

def tearDownModule():
    pass

class TestFunctions(unittest.TestCase):
    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_junction_inclusion_ratio(self):
        sample_1 = "1"
        sample_2 = "2"
        expected_output=["jir_score\tA raw count\tB raw count\t%s\n" % sample_records["header"]]
        #sample_records[sample_1]="filler"
        #sample_records[sample_2]="filler"
        expected_output.append("0.5\t1\t2\t%s\n" % (sample_records[sample_1]))
        expected_output.append("0.333333333333\t2\t3\t%s\n" % (sample_records[sample_2]))
        output = query_snaptron.junction_inclusion_ratio(args, {'samples':{sample_1:{"A":1,"B":2},sample_2:{"A":2,"B":3}}}, ["A","B"], sample_records)
        expected_outputstr = "".join(expected_output)
        outputstr = "".join(output)
        self.assertEqual(outputstr, expected_outputstr)

if __name__ == '__main__':
    unittest.main()
