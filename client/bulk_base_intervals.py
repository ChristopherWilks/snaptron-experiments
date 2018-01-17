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
import os
import argparse
import query_snaptron
import clsnapconf
from SnaptronIteratorBulk import SnaptronIteratorBulk

class IntervalProcessor():
    def _reset(self,group=None,chrm=None,start=-1,end=-1):
        self.group = group
        self.group_records = []
        self.group_chrm = chrm
        self.group_start = start
        self.group_end = end

    def __init__(self):
        self.delim = '\t'
        self.col_offset = 4
        self._reset()
        self.first = True

    #vitual base functions to be implemented by child classes
    def _write_header(self):
        pass

    def _write_records(self,new_group):
        pass

    def process(self,record):
        if self.first:
            self._write_header(record)
            self.first = False
            return
        #assume records are continguously ordered by group and then by chrm/start
        fields = record.split(self.delim)
        #empty lines
        if len(fields) <= 1:
            return
        (group, chrm, start, end) = fields[:self.col_offset]
        if group != self.group:
            if self.group is not None:
                self._write_records(group)
            self._reset(group,chrm,start,end)
        self.group_end = end
        self.group_records.append(fields)
    
    def finish(self):
        if self.group is not None:
            self._write_records(None)
            self.group = None
        

class GeneExonIntervalProcessor(IntervalProcessor):

    def _reset(self,group=None,chrm=None,start=-1,end=-1):
        IntervalProcessor._reset(self,group,chrm,start,end)
        self.prev_gene_name = self._gene_name_from_group(group)
        if self.gene_start == -1:
            self.gene_start = start
        self.gene_end = end

    def __init__(self,gene_output_file,exon_output_file):
        self.gout = gene_output_file
        self.eout = exon_output_file
        self.gene_start = -1
        self.gene_sums = []

        IntervalProcessor.__init__(self)

    def _gene_name_from_group(self,group):
        if group is None:
            return None
        return group.split(':')[0]

    def _summarize_by_exon(self):
        sums = [0.0] * len(self.group_records[0][self.col_offset:])
        for counts in self.group_records:
            sums = [sums[i]+float(count) for (i,count) in enumerate(counts[self.col_offset:])]
        return sums
    
    def _write_header(self,line):
        self.gout.write(line+"\n")
        self.eout.write(line+"\n")

    def _write_records(self,new_group):
        #assume each basic group is an exon (gene_name:chrm:exon_start:exon_end)
        #handle previous exon (doesn't include current exon's counts)
        exon_sums = self._summarize_by_exon()
        
        gene_name = self._gene_name_from_group(new_group)
        #first time through the gene sums might not have been initialized
        if len(self.gene_sums) == 0:
            self.gene_sums = [0] * len(exon_sums)
        #now check to see if we've switched genes with the current exon
        #either the first gene or we've switched genes
        if gene_name != self.prev_gene_name:
            #once more to catch the last exon
            gene_sums = [str(self.gene_sums[i]+count) for (i,count) in enumerate(exon_sums)]
            #self.gene_sums = [str(x) for x in self.gene_sums]
            self.gout.write("%s\t%s\t%s\t%s\t%s\n" % (self.prev_gene_name,self.group_chrm,self.gene_start,self.gene_end,self.delim.join(gene_sums)))
            self.gene_start = -1
            self.gene_sums = []
        else:
            #update gene sums upto and including previous exon's counts
            self.gene_sums = [self.gene_sums[i]+count for (i,count) in enumerate(exon_sums)]
        
        #now print exon
        exon_sums = [str(x) for x in exon_sums]
        self.eout.write("%s\t%s\t%s\t%s\t%s\n" % (self.group,self.group_chrm,self.group_start,self.group_end,self.delim.join(exon_sums)))
        
    
def main(args):
    (query_params_per_group, groups, endpoint) = query_snaptron.parse_query_params(args)
    gout = None
    eout = None
    if args.bulk_query_gzip:
        gout = gzip.open(args.bulk_query_file + ".snapout.genes.tsv.gz", "wb")
        eout = gzip.open(args.bulk_query_file + ".snapout.exons.tsv.gz", "wb")
    else:
        gout = open(args.bulk_query_file + ".snapout.genes.tsv", "wb")
        eout = open(args.bulk_query_file + ".snapout.exons.tsv", "wb")
    processor = GeneExonIntervalProcessor(gout,eout)
    for i in xrange(0, len(query_params_per_group), clsnapconf.BULK_LIMIT):
        sIT = SnaptronIteratorBulk(query_params_per_group[i:i+clsnapconf.BULK_LIMIT], args.datasrc, endpoint, None, processor=processor)
    processor.finish()
    gout.close()
    eout.close()

if __name__ == '__main__':
    #we can always add to these
    parser = query_snaptron.create_parser()
    args = parser.parse_args()
    if args.region is None and args.metadata is None and args.query_file is None and args.bulk_query_file is None:
        sys.stderr.write("Error: no region-related arguments passed in, exiting\n")
        parser.print_help()
        sys.exit(-1)
    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
    main(args)
