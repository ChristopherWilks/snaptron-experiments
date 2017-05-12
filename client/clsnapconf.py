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

DS_SRAV2='srav2'
DS_SRAV1='srav1'
DS_GTEX='gtex'
DS_TCGA='tcga'
SERVICE_URL='http://snaptron.cs.jhu.edu'
BUFFER_SIZE_BYTES=1024*1024
CACHE_SAMPLE_METADTA=True
#TODO use python tmp
TMPDIR='../snaptron_tmp'
FIELD_ARGS={'region':['"chr#:start-end"',str,None,'either a simple genomic region (e.g. chr1:1-1000) or a gene fusion pair (e.g. EML4-ALK)'],
      'filters':['"coverage_sum>=5&annotated=1"',str,None,'one or more junction specific thresholds/filters'], 
      'metadata':['"design_description: cortex"',str,None,'one or more sample specific filters passed to Lucene (using the Lucene query langage)'],
      'samples':['"0,1,5,200"',str,None,'one or more comma delimited sample IDs filter junctions by'],
      'contains':['1',str,None,'search for only junctions which have both their genomic coordinates within (inclusive) the specified region'],
      'exact':['1',str,None,'search for only junctions which have both their genomic coordinates matching exactly the specified region\'s start/end coordinates'],
      'either':['2',str,None,'search for only junctions which have one or the other of their genomic coordinates within (inclusive) of the specified region, 1=left coordinte;2=right coordinate']}
SAMPLE_IDS_COL=12
INTRON_ID_COL=1
INTERVAL_CHR_COL=2
INTERVAL_START_COL=3
INTERVAL_END_COL=4
MIN_EXON_SIZE = 1
LEFT_ANNOT_COL=10
RIGHT_ANNOT_COL=11
FULL_ANNOT_COL=7
