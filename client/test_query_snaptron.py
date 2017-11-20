#!/usr/bin/env python
import sys
import unittest
from argparse import Namespace
import clsnapfunc as snf
import test_data

def setUpModule():
    pass

def tearDownModule():
    pass


class TestSampleCountingPerGroup(unittest.TestCase):
    '''Test the count_samples_per_group(...) method
    which underlies many of the client high-level queries/functions
    (TS, SSC, JIR).  PSI and INTERSECTION require no special 
    handling outside of what the above three do.'''

    #we'll also need to consider the 'either' modifier as another
    #source of potential bugs

    def setUp(self):
        pass
    
    def general_test_first_and_full_processing(self, query_id, raw_file, args, group, results, first_results, full_results):
        with open(raw_file, 'r') as fin:
            snf.count_samples_per_group(args, results, fin.readline().rstrip(), group, out_fh=None)
            snf.count_samples_per_group(args, results, fin.readline().rstrip(), group, out_fh=None)
           
            #first record processed (after header) test
            self.assertEqual(first_results, results) 
        
            #keep going to the end of the file
            for line in fin:
                snf.count_samples_per_group(args, results, line.rstrip(), group, out_fh=None)

            #now check the full results data structure
            self.assertEqual(full_results, results) 
    
    def test_TS_count_samples_per_group(self):
        args = Namespace(function='ts')
        
        group = 'KCNIP4'

        results = {'samples':{},'queries':[],'exons':{'start':{},'end':{}}}
        results['groups_seen']={group:1}
        results['shared']={}
        results['annotated']={group:{}}
        results['annotated'][group][results['groups_seen'][group]]=0
        results['annotations']={group:{}}
        
        #2 tests per modifier value (2 or 1)
        for either in (2,1):
            results['either'] = either
            first_results = test_data.ts_results_first_time_through[either]
            full_results = test_data.ts_results_after_first_basic_query[either]
            raw_file = './tests/ts.raw%d' % either
            self.general_test_first_and_full_processing(either, raw_file, args, group, results, first_results, full_results)


    def test_JIR_count_samples_per_group(self):
        args = Namespace(function='jir')
        results = {'samples':{},'queries':[],'exons':{'start':{},'end':{}},'either':0}

        for (idx,group) in enumerate(list(['A_NormalTSS','B_AltTSS'])):
            first_results = test_data.jir_results_first_time_through[group]
            full_results = test_data.jir_results_after_first_basic_query[group]
            raw_file = './tests/jir.raw%d' % (int(idx)+1)
            self.general_test_first_and_full_processing(group, raw_file, args, group, results, first_results, full_results)
            

    #SSC uses the core parts of the count_samples_per_group(...) method
    def test_SSC_count_samples_per_group(self):
        args = Namespace(function='shared')
        
        group = 'C1orf100 non_validated'

        results = {'samples':{},'queries':[],'exons':{'start':{},'end':{}}}
        results['groups_seen']={group:1}
        results['shared']={}
        results['annotated']={group:{}}
        results['annotated'][group][results['groups_seen'][group]]=0
        results['annotations']={group:{}}
        
        #2 tests per modifier value (1 or 2)
        for either in (2,1):
            results['either'] = either
            first_results = test_data.C1orf100_results_first_time_through[either]
            full_results = test_data.C1orf100_results_after_first_basic_query[either]
            raw_file = './tests/C1orf100.raw%d' % either
            self.general_test_first_and_full_processing(either, raw_file, args, group, results, first_results, full_results)


class TestJIR(unittest.TestCase):
    
    def setUp(self):
        pass

    #using TCGA as the datasrc
    def test_junction_inclusion_ratio(self):
        args = Namespace(function='jir',noheader=False,limit=-1,datasrc='tcga')
        groups = ['A_NormalTSS', 'B_AltTSS'] 
        group_list = set()
        map(lambda x: group_list.add(x), groups)
        group_list = sorted(group_list)

        output = snf.junction_inclusion_ratio(args,test_data.jir_full_results,group_list,None)
        with open("./tests/jir_tcga_stats_output.tsv","rb") as fin:
            lines = fin.read()
            lines = lines.split('\n')
            lines = lines[:len(lines)-1]
            lines = map(lambda x: x+"\n", lines)
            self.assertEqual(lines, output)


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
