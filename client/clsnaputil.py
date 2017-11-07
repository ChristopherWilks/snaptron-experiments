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
import gzip

import clsnapconf

#assumes you already have the AUCs
#pulled out using:
#wiggletools print non_unique_base_coverage.bw.auc AUC non_unique_base_coverage.bw

NORMAL_RECOUNT_TARGET = 40 * 1000000

def normalize_coverage(args, record):
    #TODO: support different normalization methods, one for gene ('recount')
    #and another for junctions (use total junction coverage?)
    auc_col = clsnapconf.AUC_COL_MAP[args.datasrc]
    fields = record.rstrip().split('\t')
    if fields[1] == 'snaptron_id':
        return record
    fields[clsnapconf.SAMPLE_IDS_COL] = ','+",".join([x.split(':')[0]+":"+str((NORMAL_RECOUNT_TARGET * float(x.split(':')[1]))/float(args.sample_records_split[x.split(':')[0]][auc_col])) for x in fields[clsnapconf.SAMPLE_IDS_COL].split(',') if x != '' and x.split(':')[0] in args.sample_records_split])
    return "\t".join(fields)
