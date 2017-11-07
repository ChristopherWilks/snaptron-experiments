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
import math

import clsnapconf

#assumes you already have the AUCs
#pulled out using:
#wiggletools print non_unique_base_coverage.bw.auc AUC non_unique_base_coverage.bw

#recount's gene normalization scaling factor
NORMAL_GENE_TARGET = 40 * 1000000
#our junction scaling factor
NORMAL_JUNCTION_TARGET = 1000 * 1000000

def round_like_R(num, ndigits=0):
    '''Attempt to do IEC 60559 rounding half-way cases to nearest even (what R uses) to be equivalent to recount''' 
    p = 10**(max(0,ndigits-1))
    absx = math.fabs(num*p)
    y = math.floor(absx)
    diff = absx - y
    if diff > 0.5 or (diff == 0.5 and y % 2 != 0):
        return math.copysign((y / p) + 1.0, num)
    return math.copysign(y / p, num)


def normalize_coverage(args, record, scaling_factor=NORMAL_GENE_TARGET):
    auc_col = clsnapconf.AUC_COL_MAP[args.datasrc]
    fields = record.rstrip().split('\t')
    if fields[1] == 'snaptron_id':
        return record
    fields[clsnapconf.SAMPLE_IDS_COL] = ','+",".join([y for y in [x.split(':')[0]+":"+str(int(round_like_R((scaling_factor * float(x.split(':')[1]))/float(args.sample_records_split[x.split(':')[0]][auc_col])))) for x in fields[clsnapconf.SAMPLE_IDS_COL].split(',') if x != '' and x.split(':')[0] in args.sample_records_split] if y.split(':')[1] != "0"])
    #fields[clsnapcopnf.SAMPLE_SUM_COL] = sum([int(x.split(':')[1]) for x in fields[clsnapconf.SAMPLE_IDS_COL].split(',')])
    return "\t".join(fields)
