#!/usr/bin/env python2.7
import sys
import os
import urllib
import urllib2
import argparse
import gzip
import csv
import re

import clsnapconf
from SnaptronIteratorHTTP import SnaptronIteratorHTTP
from SnaptronIteratorLocal import SnaptronIteratorLocal

GTEX_TISSUE_COL=65

fmap = {'thresholds':'rfilter','filters':'sfilter','region':'regions'}
def parse_query_argument(record, fieldnames, groups, header=True):
    endpoint = 'snaptron'
    query=[]
    for field in fieldnames:
        if len(record[field]) > 0:
            if field == 'thresholds' or field == 'filters':
                predicates = re.sub("=",":",record[field])
                predicates = predicates.split('&')
                query.append("&".join(["%s=%s" % (fmap[field],x) for x in predicates]))
            elif field == 'group':
                groups.append(record[field])
            else:
                mapped_field = field
                if field in fmap:
                    mapped_field = fmap[field]
                query.append("%s=%s" % (mapped_field,record[field]))
    if not header:
        query.append("header=0")
    return (query,endpoint)


def parse_command_line_args(args):
    fieldnames = []
    endpoint = 'snaptron'
    for field in clsnapconf.FIELD_ARGS.keys():
        if field in vars(args) and vars(args)[field] is not None:
            fieldnames.append(field)
    groups = []
    (query,endpoint) = parse_query_argument(vars(args), fieldnames, groups, header=args.function is not None or not args.noheader)
    return (["&".join(query)], groups, endpoint)


def parse_query_params(args):
    if args.query_file is None:
        return parse_command_line_args(args)
    endpoint = 'snaptron'
    queries = []
    groups = []
    with open(args.query_file,"r") as cfin:
        creader = csv.DictReader(cfin,dialect=csv.excel_tab)
        get_header = args.function is not None or not args.noheader
        for (i,record) in enumerate(creader):
            (query, endpoint) = parse_query_argument(record, creader.fieldnames, groups, header=get_header)
            queries.append("&".join(query))
            if args.function is None:
                get_header = False
    #assume the endpoint will be the same for all lines in the file
    return (queries,groups,endpoint)

def calc_jir(a, b):
    numer = b - a
    denom = a + b + 1
    return numer/float(denom)


def junction_inclusion_ratio(args, results, group_list, sample_records):
    sample_stats = results['samples']
    group_a = group_list[0]
    group_b = group_list[1]
    sample_scores = {}
    for sample in sample_stats:
        if group_a not in sample_stats[sample]:
            sample_stats[sample][group_a]=0
        if group_b not in sample_stats[sample]:
            sample_stats[sample][group_b]=0
        sample_scores[sample] = calc_jir(sample_stats[sample][group_a], sample_stats[sample][group_b])
    missing_sample_ids = set()
    counter = 0
    if not args.noheader:
        sheader = sample_records["header"]
        sys.stdout.write("jir_score\t%s raw count\t%s raw count\t%s\n" % (group_a,group_b,sheader))
    for sample in sorted(sample_scores.keys(),key=sample_scores.__getitem__,reverse=True):
        counter += 1
        if args.limit > -1 and counter > args.limit:
            break
        score = sample_scores[sample]
        if sample not in sample_records:
            missing_sample_ids.add(sample)
            continue
        sample_record = sample_records[sample]
        sys.stdout.write("%s\t%d\t%d\t%s\n" % (str(score),sample_stats[sample][group_a],sample_stats[sample][group_b],sample_record))

def track_exons(args, results, record, group, out_fh=None):
    exons = results['exons']
    if out_fh is not None:
        out_fh.write(record+"\n")
    if 'snaptron_id' in record:
        return
    fields = record.split('\t')
    snid = fields[clsnapconf.INTRON_ID_COL]
    for (type_,col) in {'start':clsnapconf.INTERVAL_START_COL,'end':clsnapconf.INTERVAL_END_COL}.iteritems():
        coord = int(fields[col])
        if coord not in exons:
            exons[coord]={type_:set()}
        if type_ not in exons[coord]:
            exons[coord][type_]=set()
        exons[coord][type_].add(snid)

def filter_exons(args, results, group_list, sample_records):
    #only used if filtering by length range
    (rlen1,rlen2) = (None,None)
    if args.exon_length is not None:
        (rlen1,rlen2) = map(lambda x: int(x), args.exon_length.split('-'))
    #filter the joint list of both intron start and end coordinates
    coords = sorted(results['exons'].keys())
    end = None
    end_ids = None
    sys.stdout.write("Type\tLeft End Snaptron IDs\tRight End Snaptron IDs\tstart\tend\tlength\n")
    for (i,coord) in enumerate(coords):
        fields = results['exons'][coord]
        #end here is intron end (exon start)
        if 'end' in fields:
            end = coord
            end_ids = fields['end']
            #look for joining exon ends after this entry
            for coord2 in coords[i+1:]:
                fields2 = results['exons'][coord2]
                #start here is intron start (exon end)
                if 'start' in fields2:
                    start_ids = fields2['start']
                    dist = coord2 - 1 - end
                    if dist >= clsnapconf.MIN_EXON_SIZE and (rlen1 is None or (dist >= int(rlen1) and dist <= int(rlen2))):
                        sys.stdout.write("exon\t%s\t%s\t%d\t%d\t%d\n" % (",".join(end_ids),",".join(start_ids),end+1,coord2-1,dist))

TISSUE_SPECIFICITY_FUNC='ts'
SHARED_SAMPLE_COUNT_FUNC='shared'
JIR_FUNC='jir'
TRACK_EXONS_FUNC='exon'
FUNCTION_TO_TYPE={TRACK_EXONS_FUNC:'not-shared', JIR_FUNC:'not-shared', None:None, TISSUE_SPECIFICITY_FUNC:'shared',SHARED_SAMPLE_COUNT_FUNC:'shared'}
def count_samples_per_group(args, results, record, group, out_fh=None):
    sample_stats = results['samples']
    if out_fh is not None:
        out_fh.write(record+"\n")
    if 'snaptron_id' in record:
        return
    fields = record.split('\t')
    #determine annotation status of the specific splice site from "either" modifier
    annot_col = clsnapconf.LEFT_ANNOT_COL if results['either']==1 else clsnapconf.RIGHT_ANNOT_COL
    if results['either'] == 0:
        annot_col = clsnapconf.FULL_ANNOT_COL
    annotated = 1 if fields[annot_col] != "0" else 0
    if FUNCTION_TO_TYPE[args.function] == 'shared':
        results['annotated'][group][results['groups_seen'][group]]+=annotated
    #get sample list 
    samples = fields[clsnapconf.SAMPLE_IDS_COL].split(',')
    sample_covs = fields[clsnapconf.SAMPLE_IDS_COL+1].split(',')
    start_value = 0
    if args.function == TISSUE_SPECIFICITY_FUNC:
        start_value = sys.maxint
    #track samples shared across the flanking junctions here
    if 'shared' in results and group not in results['shared']:
        results['shared'][group] = set()
    for (i,sample_id) in enumerate(samples):
        #this can happen with GTEx
        if int(sample_covs[i]) == 0:
            continue
        #if we're doing tissue spec. then make sure we get shared samples, otherwise skip
        if FUNCTION_TO_TYPE[args.function] == 'shared' and results['groups_seen'][group] > 1:
            #haven't seen this sample before, so must not be shared
            if sample_id not in sample_stats or group not in sample_stats[sample_id]:
                continue
            else:
                results['shared'][group].add(sample_id)
        if sample_id not in sample_stats:
            sample_stats[sample_id]={}
        if group not in sample_stats[sample_id]:
            sample_stats[sample_id][group]=start_value
        if FUNCTION_TO_TYPE[args.function] != TISSUE_SPECIFICITY_FUNC:
            sample_stats[sample_id][group]+=int(sample_covs[i])
        else:
            #initially we used coverage of just those shared samples, but now we just do present or not
            #sample_stats[sample_id][group]=min(sample_stats[sample_id][group],int(sample_covs[i]))
            sample_stats[sample_id][group]=1

def tissue_specificity(args, results, group_list, sample_records):
    sample_stats = results['samples']
    sys.stdout.write("group\tsample_id\tshared\ttissue\n")
    for group in group_list:
        if group not in results['shared'] or len(results['shared'][group]) == 0:
            sys.stderr.write("No shared samples between splice junctions for %s\n" % (group))
        #for sample_id in results['shared'][group]:
        for sample_id in sample_records.keys():
            present = 0
            #if sample_id in sample_stats and group in sample_stats[sample_id]:
            if sample_id in results['shared'][group]:
                present = sample_stats[sample_id][group]
            sfields = sample_records[sample_id].split("\t")
            tissue = sfields[GTEX_TISSUE_COL]
            sys.stdout.write("%s\t%s\t%d\t%s\n" % (group, sample_id, present, tissue))

def report_shared_sample_counts(args, results, group_list, sample_records):
    sys.stdout.write("group\tshared_sample_counts\n")
    shared_group_count = 0
    total_fully_annotated_count = 0
    for group in group_list:
        #now see if this is a shared group or not
        if group not in results['shared'] or len(results['shared'][group]) == 0:
            sys.stderr.write("No shared samples between splice junctions for %s, skipping\n" % (group))
            #sys.stdout.write("%s\t0\n" % (group))
            continue
        count = len(results['shared'][group])
        shared_group_count+=1
        sys.stdout.write("%s\t%d\n" % (group, count))
        #determine how many fully annotated groups there are
        annot_count_across_iters = 0
        for count_per_iter in results['annotated'][group].values():
            if count_per_iter > 0:
                annot_count_across_iters+=1
        if annot_count_across_iters == results['groups_seen'][group]:
            total_fully_annotated_count+=1
    sys.stderr.write("total # of groups with shared samples:\t%d\n" % (shared_group_count))
    sys.stderr.write("total # of groups with fully annotated splices:\t%d\n" % (total_fully_annotated_count))

def download_sample_metadata(args):
    sample_records = {}
    cache_file = os.path.join(args.tmpdir,"snaptron_sample_metadata_cache.%s.tsv.gz" % args.datasrc)
    gfout = None
    if clsnapconf.CACHE_SAMPLE_METADTA:
        if os.path.exists(cache_file):
            with gzip.open(cache_file,"r") as gfin:
                for (i,line) in enumerate(gfin):
                    line = line.rstrip()
                    fields = line.split('\t')
                    sample_records[fields[0]]=line
                    if i == 0:
                        sample_records["header"]=line
            if '' in sample_records:
                del sample_records['']
            return sample_records
        else:
            gfout = gzip.open(cache_file+".tmp","w")
    response = urllib2.urlopen("%s/%s/samples?all=1" % (clsnapconf.SERVICE_URL,args.datasrc))
    all_records = response.read()
    all_records = all_records.split('\n')
    for (i,line) in enumerate(all_records):
        fields = line.split('\t')
        sample_records[fields[0]]=line
        if i == 0:
            #remove lucene index type chars from header
            line = re.sub('_[itsf]\t','\t',line)
            line = re.sub('_[itsf]$','',line)
            sample_records["header"]=line
        if gfout is not None:
            gfout.write("%s\n" % (line))
    if gfout is not None:
        gfout.close()
        os.rename(cache_file+".tmp", cache_file)
    if '' in sample_records:
        del sample_records['']
    return sample_records

def process_group(args, group_idx, groups, group_fhs, results):
    group = None
    group_fh = None
    if group_idx < len(groups):
        group = groups[group_idx]
        if not args.noraw:
            if group in group_fhs:
                group_fh = group_fhs[group]
            else:
                group_fh = group_fhs[group] = open("%s/%s.raw" % (args.tmpdir,group),"w")
                #group_fh = group_fhs[group]
        if 'groups_seen' in results:
            #track # of items in each group
            if group not in results['groups_seen']:
                results['groups_seen'][group]=0
            if group not in results['annotated']:
                results['annotated'][group]={}
            results['groups_seen'][group]+=1
            results['annotated'][group][results['groups_seen'][group]]=0
    return (group, group_fh)


iterator_map = {True:SnaptronIteratorLocal, False:SnaptronIteratorHTTP}
either_patt = re.compile(r'either=(\d)')
def process_queries(args, query_params_per_region, groups, endpoint, function=None, local=False):
    results = {'samples':{},'queries':[],'exons':{'start':{},'end':{}},'either':0}
    if FUNCTION_TO_TYPE[args.function] == 'shared':
        results['groups_seen']={}
        results['shared']={}
        results['annotated']={}
    first = True
    group_fhs = {}
    for (group_idx, query_param_string) in enumerate(query_params_per_region):
        m = either_patt.search(query_param_string)
        if m is not None:
            results['either'] = int(m.group(1))
        sIT = iterator_map[local](query_param_string, args.datasrc, endpoint)
        #assume we get a header in this case and don't count it against the args.limit
        (group, group_fh) = process_group(args, group_idx, groups, group_fhs, results)
        counter = -1
        if args.noheader:
            counter = 0
        for record in sIT:
            if function is not None:
                function(args, results, record, group, out_fh=group_fh)
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
        first = False
    for (group,group_fh) in group_fhs.iteritems():
        group_fh.close()
    return results


compute_functions={JIR_FUNC:(count_samples_per_group,junction_inclusion_ratio),TRACK_EXONS_FUNC:(track_exons,filter_exons),TISSUE_SPECIFICITY_FUNC:(count_samples_per_group,tissue_specificity),SHARED_SAMPLE_COUNT_FUNC:(count_samples_per_group,report_shared_sample_counts),None:(None,None)}
def main(args):
    #parse original set of queries
    (query_params_per_region, groups, endpoint) = parse_query_params(args)
    #get original functions (if passed in)
    (count_function, summary_function) = compute_functions[args.function]
    #process original queries
    results = process_queries(args, query_params_per_region, groups, endpoint, function=count_function, local=args.local)
    #if either the user wanted the JIR to start with on some coordinate groups, do the JIR now
    if args.function:
        sample_records = {}
        if args.function != TRACK_EXONS_FUNC:
            sample_records = download_sample_metadata(args)
        group_list = set()
        map(lambda x: group_list.add(x), groups)
        group_list = sorted(group_list)
        summary_function(args,results,group_list,sample_records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Snaptron command line client')
    for (field,settings) in clsnapconf.FIELD_ARGS.iteritems():
        parser.add_argument("--%s" % field, metavar=settings[0], type=settings[1], default=settings[2], help=settings[3])

    parser.add_argument('--query-file', metavar='/path/to/file_with_queries', type=str, default=None, help='path to a file with one query per line where a query is one or more of a region (HUGO genename or genomic interval) optionally with one or more thresholds and/or filters specified and/or contained flag turned on')

    parser.add_argument('--function', metavar='jir', type=str, default=None, help='function to compute between specified groups of junctions ranked across samples; currently only supports Junction Inclusion Ratio "jir" and exon finding "exon"')

    parser.add_argument('--tmpdir', metavar='/path/to/tmpdir', type=str, default=clsnapconf.TMPDIR, help='path to temporary storage for downloading and manipulating junction and sample records')
    
    parser.add_argument('--limit', metavar='1', type=int, default=-1, help='# of records to return, defaults to all (-1)')

    parser.add_argument('--exon-length', metavar='50-60', type=str, default=None, help='length range of exons to look for within queried region when function is set to "exon", defaults to None (print out all exons in queried region when function is "exon")')
    
    
    parser.add_argument('--local', action='store_const', const=True, default=False, help='if running Snaptron modeules locally (skipping WSI)')
    parser.add_argument('--noraw', action='store_const', const=True, default=False, help='don\'t store output of WSI in tmpdir')
    parser.add_argument('--noheader', action='store_const', const=True, default=False, help='turn off printing header in output')
    
    parser.add_argument('--datasrc', metavar='data_source_namd', type=str, default=clsnapconf.DS_SRAV2, help='data source instance of Snaptron to use, check clsnapconf.py for list')
    

    #returned format (UCSC, and/or subselection of fields) option?
    #intersection or union of intervals option?


    args = parser.parse_args()
    if args.region is None and args.thresholds is None and args.filters is None and args.query_file is None:
        sys.stderr.write("Error: no discernible arguments passed in, exiting\n")
        parser.print_help()
        sys.exit(-1)
    if args.function == TISSUE_SPECIFICITY_FUNC and 'gtex' not in args.datasrc:
        sys.stderr.write("Error: attempting to do tissue specificity (ts) on a non-GTEx Snaptron instance. Please change the SERVICE_URL setting in clsnapconf.py file to be the GTEx Snaptron instance before running this function; exiting\n")
        sys.exit(-1)
    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
    main(args)
