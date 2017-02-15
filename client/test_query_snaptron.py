#!/usr/bin/env python2.7

# Copyright 2016, Christopher Wilks <broadsword@gmail.com>
#
# This file is part of Snaptron.
#
# Snaptron is free software: you can redistribute it and/or modify
# it under the terms of the 
# Creative Commons Attribution-NonCommercial 4.0 
# International Public License ("CC BY-NC 4.0").
#
# Snaptron is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# CC BY-NC 4.0 license for more details.
#
# You should have received a copy of the CC BY-NC 4.0 license
# along with Snaptron.  If not, see 
# <https://creativecommons.org/licenses/by-nc/4.0/legalcode>.

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
        self.groupA='A'
        self.groupB='B'
        self.groupC='C'
        self.sample_1 = "1"
        self.sample_2 = "2"
        pass
    
    def tearDown(self):
        pass

    def test_junction_inclusion_ratio(self):
        expected_output=["jir_score\t%s raw count\t%s raw count\t%s\n" % (self.groupA,self.groupB,sample_records["header"])]
        #sample_records[sample_1]="filler"
        #sample_records[sample_2]="filler"
        expected_output.append("0.25\t1\t2\t%s\n" % (sample_records[self.sample_1]))
        expected_output.append("0.166666666667\t2\t3\t%s\n" % (sample_records[self.sample_2]))
        #expected_output.append("0.333333333333\t2\t3\t%s\n" % (sample_records[self.sample_2]))
        expected_outputstr = "".join(expected_output)

        output = query_snaptron.junction_inclusion_ratio(args, {'samples':{self.sample_1:{self.groupA:1,self.groupB:2},self.sample_2:{self.groupA:2,self.groupB:3}}}, [self.groupA, self.groupB], sample_records)
        outputstr = "".join(output)
        self.assertEqual(outputstr, expected_outputstr)

    def test_report_shared_sample_counts(self):
        results = {'shared':{self.groupA:set([self.sample_1]),self.groupB:([self.sample_2]),self.groupC:set([self.sample_1,self.sample_2])},
                   'annotated':{self.groupA:{1:2,2:1},self.groupB:{1:2,2:0,3:0},self.groupC:{1:0,2:0}},
                   'annotations':{self.groupA:{"sG19":[1,1],"aC19":[1,1]},self.groupB:{"sG38":[1,0]}},
                   'groups_seen':{self.groupA:2,self.groupB:3,self.groupC:2}}
        expected_output = ["group\tshared_sample_counts\n"]
        expected_output.append("%s\t%d\n" % (self.groupA, 1))
        expected_output.append("%s\t%d\t%s\n" % (self.groupA, 2, "sG19:1,1;aC19:1,1"))
        expected_output.append("%s\t%d\n" % (self.groupB, 1))
        expected_output.append("%s\t%d\n" % (self.groupC, 2))
        expected_output.append("total # of groups with shared samples:\t3\n")
        expected_output.append("total # of groups with fully annotated splices:\t1\n" )
        expected_outputstr = "".join(expected_output)
        output = query_snaptron.report_shared_sample_counts(args, results, [self.groupA, self.groupB, self.groupC], sample_records)
        outputstr = "".join(output)
        self.assertEqual(outputstr, expected_outputstr)

if __name__ == '__main__':
    unittest.main()
