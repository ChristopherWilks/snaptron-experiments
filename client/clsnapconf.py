#!/usr/bin/env python2.7

# This file is part of Snaptron.
#
# Snaptron is free software: you can redistribute it and/or modify
# it under the terms of the 
#
#    The MIT License
#
#    Copyright (c) 2016-  by Christopher Wilks <broadsword@gmail.com> 
#                         and Ben Langmead <langmea@cs.jhu.edu>
#
#    Permission is hereby granted, free of charge, to any person obtaining
#    a copy of this software and associated documentation files (the
#    "Software"), to deal in the Software without restriction, including
#    without limitation the rights to use, copy, modify, merge, publish,
#    distribute, sublicense, and/or sell copies of the Software, and to
#    permit persons to whom the Software is furnished to do so, subject to
#    the following conditions:
#
#    The above copyright notice and this permission notice shall be
#    included in all copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#    BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#    ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#    SOFTWARE.

DS_SRAV2='srav2'
DS_SRAV1='srav1'
DS_GTEX='gtex'
DS_TCGA='tcga'
DS_SUPERMOUSE='supermouse'
DS_MESA='mesa'
DS_ENCODE='encode1159'
DS_CT_H_S='ct_h_s'
DS_CT_H_B='ct_h_b'
DS_CT_M_S='ct_m_s'
DS_CT_M_B='ct_m_b'
BASE_SAMPLE_ID_OFFSETS={DS_TCGA:59761,DS_GTEX:50099,DS_SRAV2:0,DS_SRAV1:0,DS_SUPERMOUSE:0,DS_MESA:0,DS_ENCODE:0}
SERVICE_URL='https://snaptron.cs.jhu.edu'
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
SAMPLE_COUNT_COL=13
SAMPLE_SUM_COL=14
SAMPLE_AVG_COL=15
SAMPLE_MED_COL=16
INTRON_ID_COL=1
INTERVAL_CHR_COL=2
INTERVAL_START_COL=3
INTERVAL_END_COL=4
STRAND_COL=6
MIN_EXON_SIZE = 1
LEFT_ANNOT_COL=10
RIGHT_ANNOT_COL=11
FULL_ANNOT_COL=7
EXON_COUNT_COL=9
MIN_COUNT=20
BULK_LIMIT=50
BULK_QUERY_DELIMITER='|||'
ID_LIMIT=100
JX_ENDPOINT="snaptron"
SAMPLE_ENDPOINT="sample"
ANNOTATION_ENDPOINT="annotation"
GENES_ENDPOINT="genes"
EXONS_ENDPOINT="exons"
BASES_ENDPOINT="bases"
HTTP_ENDPOINTS={JX_ENDPOINT:JX_ENDPOINT,SAMPLE_ENDPOINT:'samples',ANNOTATION_ENDPOINT:'annotations',GENES_ENDPOINT:GENES_ENDPOINT,EXONS_ENDPOINT:EXONS_ENDPOINT,BASES_ENDPOINT:BASES_ENDPOINT}
LOCAL_ENDPOINTS={JX_ENDPOINT:'snaptron.py',SAMPLE_ENDPOINT:'snample.py',ANNOTATION_ENDPOINT:'snannotation.py'}
#normalization related config parameters
#recount's gene normalization scaling factor
NORMAL_RECOUNT_GENE_TARGET = 40 * 1000000
#our junction scaling factor: 40 million 100-base reads
NORMAL_RECOUNT_JX_TARGET = NORMAL_RECOUNT_GENE_TARGET * 100
#based on average junction coverage across 43332 SRAv2 samples, rounded
NORMAL_JXCOV_TARGET = 3953678
RECOUNT_NORM='recount'
JX_NORM='jxcov'
#map to the various scaling factors we can use depending on 1) normalization type and 2) endpoint (genes, junctions)
RECOUNT_FACTOR_MAP={GENES_ENDPOINT:NORMAL_RECOUNT_GENE_TARGET,JX_ENDPOINT:NORMAL_RECOUNT_JX_TARGET,EXONS_ENDPOINT:NORMAL_RECOUNT_GENE_TARGET}
JX_FACTOR_MAP={JX_ENDPOINT:NORMAL_JXCOV_TARGET}
NORM_FACTOR_MAP={RECOUNT_NORM:RECOUNT_FACTOR_MAP,JX_NORM:JX_FACTOR_MAP}
#map to the various normalization dividers depending on 1) normalization type and 2) datasource
AUC_COL_MAP={DS_SRAV2:78,DS_GTEX:308,DS_TCGA:11,DS_SUPERMOUSE:55,DS_MESA:55,DS_ENCODE:37,DS_CT_H_S:51,DS_CT_H_B:51,DS_CT_M_S:51,DS_CT_M_B:51}
JXCOV_COL_MAP={DS_SRAV2:90,DS_GTEX:320,DS_TCGA:859,DS_SUPERMOUSE:53,DS_MESA:53,DS_ENCODE:39,DS_CT_H_S:49,DS_CT_H_B:49,DS_CT_M_S:49,DS_CT_M_B:49}
NORM_DIVISOR_COL_MAP={RECOUNT_NORM:AUC_COL_MAP,JX_NORM:JXCOV_COL_MAP}
