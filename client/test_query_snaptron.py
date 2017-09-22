#!/usr/bin/env python2.7
import sys
import unittest
from argparse import Namespace
import query_snaptron

def setUpModule():
    pass

def tearDownModule():
    pass


class TestSampleCountingPerGroup(unittest.TestCase):
    '''Test the count_samples_per_group(...) method
    which underlies many of the client high-level queries/functions
    (TS, SSC, JIR, PSI)'''

    #we'll also need to consider the 'either' modifier as another
    #source of potential bugs

    def setUp(self):
        pass
    
    def test_ssc_count_samples_per_group(self):
        group = 'C1orf100 non_validated'
        either = 2
        results = {'samples':{},'queries':[],'exons':{'start':{},'end':{}},'either':either}
        results['groups_seen']={group:1}
        results['shared']={}
        results['annotated']={group:{}}
        results['annotated'][group][results['groups_seen'][group]]=0
        results['annotations']={group:{}}

        args = Namespace(function='shared')

        with open('../tests/C1orf100.raw','r') as fin:
            query_snaptron.count_samples_per_group(args, results, fin.readline().rstrip(), group, out_fh=None)
            query_snaptron.count_samples_per_group(args, results, fin.readline().rstrip(), group, out_fh=None)

        first_time_through_results = {'annotated': {'C1orf100 non_validated': {1: 0}}, 'exons': {'start': {}, 'end': {}}, 'groups_seen': {'C1orf100 non_validated': 1}, 'either': 2, 'samples': {'1588': {'C1orf100 non_validated': 3}, '41135': {'C1orf100 non_validated': 2}, '16768': {'C1orf100 non_validated': 1}, '18048': {'C1orf100 non_validated': 2}, '13573': {'C1orf100 non_validated': 1}, '39148': {'C1orf100 non_validated': 1}, '26599': {'C1orf100 non_validated': 1}, '7538': {'C1orf100 non_validated': 4}, '3520': {'C1orf100 non_validated': 1}, '30111': {'C1orf100 non_validated': 1}}, 'queries': [], 'shared': {'C1orf100 non_validated': set([])}, 'annotations': {'C1orf100 non_validated': {}}}
        self.assertEqual(first_time_through_results, results) 


class TestJIR(unittest.TestCase):
    
    def setUp(self):
        pass

class TestSSC(unittest.TestCase):
    
    def setUp(self):
        pass

class TestTS(unittest.TestCase):
    
    def setUp(self):
        pass


#self.assertEqual(iids, EXPECTED_IIDS[IQs[i]+str(RQs_flat[r])+str(RQs_flat[r+1])])
class TestBasicQueries(unittest.TestCase):
    def setUp(self):
        pass

if __name__ == '__main__':
    unittest.main()
