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

import sys
import os
import urllib
import urllib2
import argparse
import gzip
import csv
import re
from operator import itemgetter

import clsnapconf
import clsnaputil
import clsnapfunc as snf
from SnaptronIteratorHTTP import SnaptronIteratorHTTP
from SnaptronIteratorLocal import SnaptronIteratorLocal
from SnaptronIteratorBulk import SnaptronIteratorBulk

csv.field_size_limit(sys.maxsize)

compute_functions={snf.MATES_FUNC:(snf.sum_sample_coverage,snf.report_splice_mates),clsnaputil.PSI_FUNC:(snf.count_samples_per_group,snf.percent_spliced_in),snf.JIR_FUNC:(snf.count_samples_per_group,snf.junction_inclusion_ratio),snf.TRACK_EXONS_FUNC:(snf.track_exons,snf.filter_exons),snf.TISSUE_SPECIFICITY_FUNC:(snf.count_samples_per_group,snf.tissue_specificity),snf.SHARED_SAMPLE_COUNT_FUNC:(snf.count_samples_per_group,snf.report_shared_sample_counts),snf.INTERSECTION_FUNC:(None,None),snf.GROUP_COV_FUNC:(snf.count_samples_per_group,snf.group_coverage),None:(None,None)}


def parse_query_params(args, dont_join_query_params=False):
    '''Determines whether the query was passed in via the command line
    or a file and handles the arguments appropriately'''

    if args.query_file is None and args.bulk_query_file is None:
        return clsnaputil.parse_command_line_args(args)
    endpoint = args.endpoint
    queries = []
    groups = []
    groups_seen = {}
    datasrcs = []
    infile = args.query_file
    bulk = False
    if args.bulk_query_file is not None:
        infile = args.bulk_query_file
        bulk = True
    with open(infile,"r") as cfin:
        creader = csv.DictReader(cfin,dialect=csv.excel_tab)
        get_header = args.function is not None or not args.noheader
        for (i,record) in enumerate(creader):
            (query, endpoint, datasrc) = clsnaputil.parse_query_argument(args, record, creader.fieldnames, groups, groups_seen, inline_group=bulk, header=get_header)
            if datasrc is not None:
                datasrcs.append(datasrc)
            if dont_join_query_params:
                queries.append(query)
            else:
                queries.append("&".join(query))
            if args.function is None:
                get_header = False
    if len(datasrcs) == 0:
        datasrcs.append(args.datasrc)
    #assume the endpoint will be the same for all lines in the file
    return (queries,groups,endpoint,datasrcs)


def merge(record1, record2):
    if record2 is None:
        return None
    (r1, r2) = (record1, record2)
    (group1, sid1, c1, s1, e1, len1, st1) = r1[:clsnapconf.STRAND_COL+1]
    (group2, sid2, c2, s2, e2, len2, st2) = r2[:clsnapconf.STRAND_COL+1]
    if group2 == group1 and c2 == c1 and s2 == s1 and e2 == e1 and st2 == st1:
        #merge the two records
        r1[1]+=":"+r2[1]
        #header
        if c1 == 'chromosome':
            return r1
        r1[clsnapconf.SAMPLE_IDS_COL] += r2[clsnapconf.SAMPLE_IDS_COL]
        r1[clsnapconf.SAMPLE_COUNT_COL] = str(int(r1[clsnapconf.SAMPLE_COUNT_COL]) + int(r2[clsnapconf.SAMPLE_COUNT_COL]))
        r1[clsnapconf.SAMPLE_SUM_COL] = str(int(r1[clsnapconf.SAMPLE_SUM_COL]) + int(r2[clsnapconf.SAMPLE_SUM_COL]))
        r1[clsnapconf.SAMPLE_AVG_COL] = str(float(r1[clsnapconf.SAMPLE_SUM_COL])/float(r1[clsnapconf.SAMPLE_COUNT_COL]))
        #recalculate median
        idx1 = int(r1[clsnapconf.SAMPLE_COUNT_COL]) / 2.0
        idx2 = idx1 + 1
        #odd
        if int(idx1) < idx1:
            idx1 = int(idx1) + 1
            idx2 = idx1
        idx1-=1
        idx2-=1
        counts = [0 if i != idx1 and i != idx2 else x for (i,x) in enumerate(sorted([int(y.split(':')[1]) for y in r1[clsnapconf.SAMPLE_IDS_COL].split(',')[1:]]))]
        r1[clsnapconf.SAMPLE_MED_COL] = str(float(counts[int(idx1)]+counts[int(idx2)])/2.0)
        r1[clsnapconf.SAMPLE_MED_COL+1] += ',' + r2[clsnapconf.SAMPLE_MED_COL+1]
        return r1
    return None

def process_combined_datasource_bulk_queries_multi(queries, datasrcs, endpoint, outfile, groups):
    '''Merge or Intersect junction records from multiple compilations'''
    records = {}
    iTs = []
    datasrc_jxs = []
    for (i,datasrc) in enumerate(datasrcs):
        iTs.append(SnaptronIteratorBulk(queries, datasrcs[i], endpoint, None))
        datasrc_jxs.append(set())
    its_done = 0
    num_its = len(iTs)
    header = True
    first_datasrc = True
    intersected_jxs = set()
    while its_done < num_its:
        for (i, sIT) in enumerate(iTs):
            if sIT is None:
                continue
            try:
                record = sIT.next()
                r1 = record.split('\t')
                (group1, sid1, c1, s1, e1, len1, st1) = r1[:clsnapconf.STRAND_COL+1]
                if c1 == 'chromosome':
                    if header:
                        outfile.write(record+'\n')
                        header = False
                    continue
                key = (group1, c1, int(s1), int(e1))
                datasrc_jxs[i].add(key)
                final_record = record
                if key in records:
                    r2 = records[key].split('\t')
                    merged = merge(r2, r1)
                    if merged is not None:
                        final_record = '\t'.join(merged)
                    records[key] = final_record
                    continue
                records[key] = final_record
            except StopIteration, si:
                its_done += 1
                iTs[i] = None
    if args.function == snf.INTERSECTION_FUNC:
        intersected_jxs = datasrc_jxs[0].intersection(*datasrc_jxs)
    #now read out all the junction rows sorted on (group, chromosome, start, end)
    keys = intersected_jxs
    if args.function != snf.INTERSECTION_FUNC:
        keys = records.keys()
    for key in sorted(keys, key=itemgetter(0,1,2,3)):
        outfile.write(records[key]+"\n")


def process_bulk_queries(args):
    (query_params_per_group, groups, endpoint, datasrcs) = parse_query_params(args)
    outfile = sys.stdout
    if not args.bulk_query_stdout:
        if args.bulk_query_gzip:
            outfile = gzip.open(args.bulk_query_file + ".snap_results.tsv.gz", "wb")
        else:
            outfile = open(args.bulk_query_file + ".snap_results.tsv", "wb")
    #we ONLY support a list of datasrcs passed in at the args level
    #i.e. ignoring an per-query datasrc passed via a query file
    if len(datasrcs) != 1:
        sys.stderr.write("additional datasources supplied via query file, ignoring and only using these: %s\n" % (args.datasrc))
    datasrcs = args.datasrc.split(',')
    #TODO: make gzip optional for output file
    for i in xrange(0, len(query_params_per_group), clsnapconf.BULK_LIMIT):
        if len(datasrcs) > 1:
            process_combined_datasource_bulk_queries_multi(query_params_per_group[i:i+clsnapconf.BULK_LIMIT], datasrcs, endpoint, outfile, groups)
        else:
            sIT = SnaptronIteratorBulk(query_params_per_group[i:i+clsnapconf.BULK_LIMIT], args.datasrc, endpoint, outfile)
    outfile.close()


def process_group(args, group_idx, groups, group_fhs, results):
    '''General method called by process_queries to handle each junction group's results for a single basic query'''
    group = None
    group_fh = None
    if group_idx < len(groups):
        group = groups[group_idx]
        if not args.noraw:
            if group in group_fhs:
                group_fh = group_fhs[group]
            else:
                group_fh = group_fhs[group] = open("%s/%s.raw" % (args.tmpdir,group),"w")
        if 'groups_seen' in results:
            #track # of items in each group
            if group not in results['groups_seen']:
                results['groups_seen'][group]=0
            if group not in results['annotated']:
                results['annotated'][group]={}
                results['annotations'][group]={}
            results['groups_seen'][group]+=1
            results['annotated'][group][results['groups_seen'][group]]=0
    return (group, group_fh)

iterator_map = {True:SnaptronIteratorLocal, False:SnaptronIteratorHTTP}
either_patt = re.compile(r'either=(\d)')
def process_queries(args, query_params_per_region, groups, endpoint, datasrcs, count_function=None, local=False):
    '''General function to process the high level queries (functions) via
    one or more basic queries while tracking the results across basic queries'''

    results = {'samples':{},'queries':[],'exons':{'start':{},'end':{}},'either':0}
    if snf.FUNCTION_TO_TYPE[args.function] == 'shared':
        results['groups_seen']={}
        results['shared']={}
        results['annotated']={}
        results['annotations']={}
    #used in intersected queries
    if snf.FUNCTION_TO_TYPE[args.function] == 'track-queries':
        results['jids']=set()
    first = True
    group_fhs = {}
    #used for intersection
    intersection_datasrc = 'srav2'
    for (group_idx, query_param_string) in enumerate(query_params_per_region):
        datasrc = datasrcs[0]
        if len(datasrcs) > 1:
            datasrc = datasrcs[group_idx]
        jids = set()
        m = either_patt.search(query_param_string)
        if m is not None:
            results['either'] = int(m.group(1))
        sIT = iterator_map[local](query_param_string, datasrc, endpoint)
        #since this will be used for intersection, doesn't matter which datasrc
        #we save as long as it's one of the ones passed in
        intersection_datasrc = datasrc
        results['siterator'] = iterator_map[local]
        #assume we get a header in this case and don't count it against the args.limit
        (group, group_fh) = process_group(args, group_idx, groups, group_fhs, results)
        counter = -1
        if args.noheader:
            counter = 0
        norm_scaling_factor = None
        norm_divisor_col = None
        if args.normalize is not None:
            norm_scaling_factor = clsnapconf.NORM_FACTOR_MAP[args.normalize][args.endpoint]
            norm_divisor_col = clsnapconf.NORM_DIVISOR_COL_MAP[args.normalize][args.datasrc]
        for record in sIT:
            #first check to see if 1) we're doing a gene expression query
            #and 2) we have a exon count constraint
            if args.exon_count > 0 and args.endpoint == clsnapconf.GENES_ENDPOINT:
                fields_ = record.split('\t')
                if int(fields_[clsnapconf.EXON_COUNT_COL]) < args.exon_count:
                    continue
            if args.normalize is not None:
                record = clsnaputil.normalize_coverage(args, record, norm_divisor_col, norm_scaling_factor)
                if record is None:
                    continue
            if count_function is not None:
                count_function(args, results, record, group, out_fh=group_fh)
            elif args.function == snf.INTERSECTION_FUNC:
                fields_ = record.split('\t')
                if fields_[clsnapconf.INTRON_ID_COL] != 'snaptron_id':
                    jids.add(int(fields_[clsnapconf.INTRON_ID_COL]))
            else:
                counter += 1
                if args.limit > -1 and counter > args.limit:
                    break
                group_label = ''
                if len(groups) > 0 and len(groups[group_idx]) > 0:
                    group_label = "%s\t" % (groups[group_idx])
                    if not args.noheader and first and counter == 0:
                        group_label = 'group\t'
                sys.stdout.write("%s%s\n" % (group_label, record))
        #keep intersecting the sets of junctions from each individual query
        if args.function == snf.INTERSECTION_FUNC:
            if first:
                results['jids'] = jids
            else:
                results['jids'] = results['jids'].intersection(jids)
        first = False
    #check if we have jids to query for (currently only in case of intersection function)
    if 'jids' in results and len(results['jids']) > 0:
        query_param_strings = clsnaputil.breakup_junction_id_query(results['jids'])
        for query_param_string in query_param_strings:
            sIT = iterator_map[local](query_param_string, datasrc, endpoint)
            for record in sIT:
                #ignore limits and group label on this output
                sys.stdout.write(record + "\n")
    for (group,group_fh) in group_fhs.iteritems():
        group_fh.close()
    return results
    

def main(args):
    sample_records = {}
    if (args.function is not None and args.function != snf.TRACK_EXONS_FUNC and args.bulk_query_file is None) or args.normalize is not None:
        (sample_records, sample_records_split) = clsnaputil.download_sample_metadata(args,split=args.normalize is not None)
        args.sample_records_split = sample_records_split
    #special handling for bulk, should attempt to refactor this in the future to
    #be more streamlined
    if args.bulk_query_file is not None:
        process_bulk_queries(args)
        return
    #parse original set of queries
    (query_params_per_region, groups, endpoint, datasrcs) = parse_query_params(args, dont_join_query_params=(args.function == clsnaputil.APSI_FUNC))
    if args.function == clsnaputil.APSI_FUNC:
        snf.process_apsi(args, query_params_per_region, groups, endpoint, datasrcs, iterator_map, local=args.local)
        return
    #get original functions (if passed in)
    (count_function, summary_function) = compute_functions[args.function]
    #process original queries
    results = process_queries(args, query_params_per_region, groups, endpoint, datasrcs, count_function=count_function, local=args.local)
    #if either the user wanted the JIR to start with on some coordinate groups, do the JIR now
    if args.function is not None and summary_function is not None:
        group_list = set()
        map(lambda x: group_list.add(x), groups)
        group_list = sorted(group_list)
        output = summary_function(args,results,group_list,sample_records)

def create_parser(disable_header=False):
    parser = argparse.ArgumentParser(description='Snaptron command line client')
    for (field,settings) in clsnapconf.FIELD_ARGS.iteritems():
        parser.add_argument("--%s" % field, metavar=settings[0], type=settings[1], default=settings[2], help=settings[3])

    parser.add_argument('--query-file', metavar='/path/to/file_with_queries', type=str, default=None, help='path to a file with one query per line where a query is one or more of a region (HUGO genename or genomic interval) optionally with one or more filters and/or metadata contraints specified and/or contained/either/within flag(s) turned on')
    
    parser.add_argument('--bulk-query-file', metavar='/path/to/file_with_queries', type=str, default=None, help='Same format as --query-file but gets run with up to 50 queries at a time (better for lots of queries).  This is the path to a file with one query per line where a query is one or more of a region (HUGO genename or genomic interval) optionally with one or more filters and/or metadata contraints specified and/or contained/either/within flag(s) turned on')
    
    parser.add_argument('--sample-group-file', metavar='/path/to/file_with_sample_groups.tsv', type=str, default=None, help='path to file with sample groupings, e.g. GTEx tissues or TCGA disease types')
    
    parser.add_argument('--bulk-query-gzip', action='store_const', const=True, default=False, help='gzip bulk query output file')

    parser.add_argument('--bulk-query-stdout', action='store_const', const=True, default=False, help='dump output of a bulk query to STDOUT instead of writing to a file')

    parser.add_argument('--function', metavar='jir', type=str, default=None, help='function to compute between specified groups of junctions ranked across samples: "jir", "psi", "apsi", "ts", "ssc", "mates", and "exon"')

    parser.add_argument('--tmpdir', metavar='/path/to/tmpdir', type=str, default=clsnapconf.TMPDIR, help='path to temporary storage for downloading and manipulating junction and sample records')
    
    parser.add_argument('--min-count', metavar='count_as_int', type=int, default=clsnapconf.MIN_COUNT, help='minimum # of reads required in total (denominator) when calculating the PSI, COVERAGE, or splice mates (--function mates, --donor, or --acceptor), default is ' + str(clsnapconf.MIN_COUNT))
    
    parser.add_argument('--min-count-jir', metavar='count_as_int', type=int, default=1, help='minimum # of total reads required across both splice groups in each sample when calculating the JIR (--function jir), default is 1')

    parser.add_argument('--summarize', action='store_const', const=True, default=False, help='report cross-sample statistics for certain functions in addition to their normal, per-sample output (e.g. PSI)')

    parser.add_argument('--limit', metavar='1', type=int, default=-1, help='# of records to return, defaults to all (-1)')

    parser.add_argument('--exon-length', metavar='50-60', type=str, default=None, help='length range of exons to look for within queried region when function is set to "exon", defaults to None (print out all exons in queried region when function is "exon")')
    parser.add_argument('--local', action='store_const', const=True, default=False, help='if running Snaptron modeules locally (skipping WSI)')
    parser.add_argument('--noraw', action='store_const', const=True, default=False, help='don\'t store output of WSI in tmpdir')
    parser.add_argument('--noheader', action='store_const', const=True, default=disable_header, help='turn off printing header in output')

    parser.add_argument('--datasrc', metavar='data_source_name', type=str, default=clsnapconf.DS_SRAV2, help='data source instance of Snaptron to use, check clsnapconf.py for list')
    parser.add_argument('--endpoint', metavar='endpoint_string', type=str, default='snaptron', help='endpoint to use ["%s"=>junctions, "%s"=>gene expression, "%s"=>exon expression, "%s"=>base expression]' % (clsnapconf.JX_ENDPOINT, clsnapconf.GENES_ENDPOINT, clsnapconf.EXONS_ENDPOINT, clsnapconf.BASES_ENDPOINT))
    parser.add_argument('--exon-count', metavar='5', type=int, default=0, help='number of exons required for any gene returned in gene expression queries, otherwise ignored')
    parser.add_argument('--normalize', metavar='normalization_string', type=str, default=None, help='Normalize? and if so what method to use [None,%s,%s]' % (clsnapconf.RECOUNT_NORM, clsnapconf.JX_NORM))
    parser.add_argument('--donor', metavar='strand_string', type=str, default=None, help='shorthand for "mates" function starting with a donor on a particular strand [+|-]; must specify a region for the splice site via --region')
    parser.add_argument('--acceptor', metavar='strand_string', type=str, default=None, help='shorthand "mates" function starting with an acceptor on a particular strand [+|-]; must specify a region for the splice site via --region')
    parser.add_argument('--event-type', metavar='event_type_string', type=str, default=None, help='only used with when looking for splice mates; [ri=retained_intron], default is None')
    return parser

splice_mates_map={'d+':'1','d-':'2','a+':'2','a-':'1'}
if __name__ == '__main__':
    parser = create_parser()
    #returned format (UCSC, and/or subselection of fields) option?
    #intersection or union of intervals option?

    args = parser.parse_args()
    if args.region is None and args.metadata is None and args.query_file is None and args.bulk_query_file is None:
        sys.stderr.write("Error: no region-related arguments passed in, exiting\n")
        parser.print_help()
        sys.exit(-1)
    if args.function == snf.TISSUE_SPECIFICITY_FUNC and 'gtex' not in args.datasrc:
        sys.stderr.write("Error: attempting to do tissue specificity (ts) on a non-GTEx Snaptron instance. Please change the SERVICE_URL setting in clsnapconf.py file to be the GTEx Snaptron instance before running this function; exiting\n")
        sys.exit(-1)
    if (args.donor or args.acceptor) and not args.region:
        sys.stderr.write("Error: no region argument passed in, but one of --donor or --acceptor was specified; exiting\n")
        parser.print_help()
        sys.exit(-1)
    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
    if args.donor:
        args.function = snf.MATES_FUNC
        strand = args.donor
        args.filters='strand=%s' % strand
        args.either=splice_mates_map['d'+strand]
    if args.acceptor:
        args.function = snf.MATES_FUNC
        strand = args.acceptor
        args.filters='strand=%s' % strand
        args.either=splice_mates_map['a'+strand]
    main(args)
