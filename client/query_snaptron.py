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
import urllib
import urllib2
import argparse
import gzip
import csv
import re

import clsnapconf
from SnaptronIteratorHTTP import SnaptronIteratorHTTP
from SnaptronIteratorLocal import SnaptronIteratorLocal
from SnaptronIteratorBulk import SnaptronIteratorBulk

GTEX_TISSUE_COL=65

TISSUE_SPECIFICITY_FUNC='ts'
SHARED_SAMPLE_COUNT_FUNC='shared'
JIR_FUNC='jir'
PSI_FUNC='psi'
TRACK_EXONS_FUNC='exon'
INTERSECTION_FUNC='intersection'

fmap = {'filters':'rfilter','metadata':'sfilter','region':'regions','samples':'sids'}
def parse_query_argument(args, record, fieldnames, groups, groups_seen, inline_group=False, header=True):
    '''Called from parse_command_line_args;
    builds the Snaptron query string from one
    or more of the separate query arguments passed in fieldnames:
    region (range), filters (rfilter), metadata (sfilter),
    and samples (sids)'''

    endpoint = 'snaptron'
    query=[None]
    fields_seen = set()
    group = None
    for field in fieldnames:
        if len(record[field]) > 0:
            fields_seen.add(field)
            if field == 'filters' or field == 'metadata':
                predicates = re.sub("=",":",record[field])
                predicates = predicates.split('&')
                query.append("&".join(["%s=%s" % (fmap[field],x) for x in predicates]))
            elif field == 'group':
                group = record[field]
                #dont want to print the header multiple times for the same group
                if group in groups_seen:
                    header = False
                    if args.function == PSI_FUNC:
                        gidx = groups_seen[group]
                        groups[gidx] = "A1_" + group
                        group = "A2_" + group
                else:
                    groups_seen[group]=len(groups)
                groups.append(group)
            else:
                mapped_field = field
                if field in fmap:
                    mapped_field = fmap[field]
                query.append("%s=%s" % (mapped_field,record[field]))
    if len(fields_seen) == 1 and "metadata" in fields_seen:
        endpoint = 'sample'
    if not header:
        query.append("header=0")
    #either we got a group or we have to shift the list over by one
    if inline_group and group is not None:
        query[0] = "group=%s" % (group)
    else:
        query = query[1:]
    return (query,endpoint)


def parse_command_line_args(args):
    '''Loop through arguments passed in on the command line and parse them'''

    fieldnames = []
    endpoint = 'snaptron'
    for field in clsnapconf.FIELD_ARGS.keys():
        if field in vars(args) and vars(args)[field] is not None:
            fieldnames.append(field)
    groups = []
    (query,endpoint) = parse_query_argument(args, vars(args), fieldnames, groups, {}, header=args.function is not None or not args.noheader)
    return (["&".join(query)], groups, endpoint)


def parse_query_params(args):
    '''Determines whether the query was passed in via the command line
    or a file and handles the arguments appropriately'''

    if args.query_file is None and args.bulk_query_file is None:
        return parse_command_line_args(args)
    endpoint = 'snaptron'
    queries = []
    groups = []
    groups_seen = {}
    infile = args.query_file
    bulk = False
    if args.bulk_query_file is not None:
        infile = args.bulk_query_file
        bulk = True
    with open(infile,"r") as cfin:
        creader = csv.DictReader(cfin,dialect=csv.excel_tab)
        get_header = args.function is not None or not args.noheader
        for (i,record) in enumerate(creader):
            (query, endpoint) = parse_query_argument(args, record, creader.fieldnames, groups, groups_seen, inline_group=bulk, header=get_header)
            queries.append("&".join(query))
            if args.function is None:
                get_header = False
    #assume the endpoint will be the same for all lines in the file
    return (queries,groups,endpoint)

def process_bulk_queries(args):
    (query_params_per_group, groups, endpoint) = parse_query_params(args)
    #TODO: make gzip optional for output file
    with open(args.bulk_query_file + ".snap_results.tsv", "w") as outfile:
        for i in xrange(0, len(query_params_per_group), clsnapconf.BULK_LIMIT):
            sIT = SnaptronIteratorBulk(query_params_per_group[i:i+clsnapconf.BULK_LIMIT], args.datasrc, endpoint, outfile)


def calc_psi(args, sample_stat, group_list):
    '''Calculates the simple PSI between 2 junction 
    variants where one variant has 2 junctions (inclusion),
    and the other has only one junction (exclusion).'''

    ss = sample_stat
    gl = group_list
    (inclusion1,inclusion2,exclusion) = (sample_stat[group_list[0]],sample_stat[group_list[1]],sample_stat[group_list[2]])
    mean_inclusion = (inclusion1 + inclusion2) / 2.0
    total = mean_inclusion + exclusion
    psi = mean_inclusion / float(total)
    if inclusion1 == 0 or inclusion2 == 0 or total < args.psi_min_reads:
        psi = -1.0
    return ([psi,total,inclusion1,inclusion2,exclusion])
                

def percent_spliced_in(args, results, group_list, sample_records):
    '''Calculates the PSI across basic queries for
    for the groups of junctions specified in group_list for samples
    from sample_records; junction group coverages are reported in results per sample'''

    sample_stats = results['samples']
    #only expect 3 groups for PSI, and order of A1,A2,B
    #where A1 and A2 are inclusion groups, B is exclusion
    sample_scores = {}
    for sample in sample_stats:
        for group in group_list:
            if group not in sample_stats[sample]:
                sample_stats[sample][group]=0
        sample_scores[sample] = calc_psi(args, sample_stats[sample], group_list)
    missing_sample_ids = set()
    counter = 0
    output = []
    if not args.noheader:
        rawstr = ""
        for group in group_list:
            rawstr += "%s raw count\t" % group
        sheader = sample_records["header"]
        outstr = "PSI\ttotal_count\t%s%s\n" % (rawstr,sheader)
        output.append(outstr)
        sys.stdout.write(outstr)
    #sort by PSI
    for sample in sorted(sample_scores.keys(),key=lambda x: sample_scores[x][:2],reverse=True):
        counter += 1
        if args.limit > -1 and counter > args.limit:
            break
        scores = sample_scores[sample]
        if sample not in sample_records:
            missing_sample_ids.add(sample)
            continue
        rawstr = "\t".join([str(x) for x in sample_scores[sample][1:]])
        sample_record = sample_records[sample]
        outstr = "%s\t%s\t%s\n" % (str(sample_scores[sample][0]),rawstr,sample_record)
        output.append(outstr)
        sys.stdout.write(outstr)
    return output


def calc_jir(a, b):
    '''Short method to do the actual junctional inclusion ratio calculation
    a and b are the junction group coverages'''

    numer = (b - a)
    denom = (a + b) + 1
    return numer/float(denom)


def junction_inclusion_ratio(args, results, group_list, sample_records):
    '''Calculates the junction inclusion ratio across basic queries for
    for the groups of junctions specified in group_list for samples
    from sample_records; junction group coverages are reported in results per sample'''

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
    output = []
    if not args.noheader:
        sheader = sample_records["header"]
        outstr = "jir_score\t%s raw count\t%s raw count\t%s\n" % (group_a,group_b,sheader)
        output.append(outstr)
        sys.stdout.write(outstr)
    for sample in sorted(sample_scores.keys(),key=sample_scores.__getitem__,reverse=True):
        counter += 1
        if args.limit > -1 and counter > args.limit:
            break
        score = sample_scores[sample]
        if sample not in sample_records:
            missing_sample_ids.add(sample)
            continue
        sample_record = sample_records[sample]
        outstr = "%s\t%d\t%d\t%s\n" % (str(score),sample_stats[sample][group_a],sample_stats[sample][group_b],sample_record)
        output.append(outstr)
        sys.stdout.write(outstr)
    return output

def intersect_junctions(args, results, record, group, out_fh=None):
    '''Intersect current snaptron_ids (junction ids) with previous queries results.
       If this is the last query in the group, do intersection but then print out results.
    '''
    fields = record.split('\t')
    snid = fields[clsnapconf.INTRON_ID_COL]
    return

def track_exons(args, results, record, group, out_fh=None):
    '''Store exons across basic queries for high level exon finding query'''

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
    '''Filter by exon length'''

    #only used if filtering by length range
    (rlen1,rlen2) = (None,None)
    if args.exon_length is not None:
        (rlen1,rlen2) = map(lambda x: int(x), args.exon_length.split('-'))
    #filter the joint list of both intron start and end coordinates
    coords = sorted(results['exons'].keys())
    end = None
    end_ids = None
    output = []
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
    return output


FUNCTION_TO_TYPE={TRACK_EXONS_FUNC:'not-shared', JIR_FUNC:'not-shared', PSI_FUNC:'not-shared', None:None, TISSUE_SPECIFICITY_FUNC:'shared',SHARED_SAMPLE_COUNT_FUNC:'shared',INTERSECTION_FUNC:'track-queries'}
def count_samples_per_group(args, results, record, group, out_fh=None):
    '''Tracks shared status of samples which appear across basic queries (junctions)
    as well as annotation status of junctions; organized by junction group'''

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
        if results['either'] > 0 and annotated == 1:
            for annot_source in fields[annot_col].split(","):
                if annot_source not in results['annotations'][group]:
                    results['annotations'][group][annot_source] = [0,0]
                results['annotations'][group][annot_source][results['either']-1] = 1
    #get sample list
    samples_ = fields[clsnapconf.SAMPLE_IDS_COL][1:].split(',')
    samples = [x.split(":")[0] for x in samples_]
    sample_covs = [x.split(":")[1] for x in samples_]
    start_value = 0
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
        if args.function != TISSUE_SPECIFICITY_FUNC:
            sample_stats[sample_id][group]+=int(sample_covs[i])
        else:
            sample_stats[sample_id][group]=1


def tissue_specificity(args, results, group_list, sample_records):
    '''Process results across basic queries to get tissue specificity values'''

    sample_stats = results['samples']
    output = []
    sys.stdout.write("group\tsample_id\tshared\ttissue\n")
    for group in group_list:
        if group not in results['shared'] or len(results['shared'][group]) == 0:
            sys.stderr.write("No shared samples between splice junctions for %s\n" % (group))
        #for sample_id in results['shared'][group]:
        for sample_id in sample_records.keys():
            ts_val = 0
            #if sample_id in sample_stats and group in sample_stats[sample_id]:
            if sample_id in results['shared'][group]:
                ts_val = sample_stats[sample_id][group]
            sfields = sample_records[sample_id].split("\t")
            tissue = sfields[GTEX_TISSUE_COL]
            sys.stdout.write("%s\t%s\t%d\t%s\n" % (group, sample_id, ts_val, tissue))
    return output

def report_shared_sample_counts(args, results, group_list, sample_records):
    '''Outputs 1) shared sample counts per junction group
    2) annotation status and sources per junction group
    3) total # of groups with shared samples across member junctions
    4) total # of groups whose member junctions are annotated'''

    output = []
    outputstr = "group\tshared_sample_counts\n"
    output.append(outputstr)
    sys.stdout.write(outputstr)
    shared_group_count = 0
    total_fully_annotated_count = 0
    annots_fh = None
    if not args.noraw:
        annots_fh = open("%s/fully_annotated_groups.%s.tsv" % (args.tmpdir,args.datasrc),"w")
    for group in group_list:
        #now see if this is a shared group or not
        if group not in results['shared'] or len(results['shared'][group]) == 0:
            outputstr = "No shared samples between splice junctions for %s, skipping\n" % (group)
            output.append(outputstr)
            sys.stderr.write(outputstr)
            continue
        count = len(results['shared'][group])
        shared_group_count+=1
        outputstr = "%s\t%d\n" % (group, count)
        output.append(outputstr)
        sys.stdout.write(outputstr)
        #determine how many fully annotated groups there are
        annot_count_across_iters = 0
        for count_per_iter in results['annotated'][group].values():
            if count_per_iter > 0:
                annot_count_across_iters+=1
        if annot_count_across_iters == results['groups_seen'][group]:
            total_fully_annotated_count+=1
            annot_sources = ";".join([x+":"+str(y[0])+","+str(y[1]) for (x,y) in results['annotations'][group].iteritems()])
            outputstr = "%s\t%d\t%s\n" % (group, annot_count_across_iters, annot_sources)
            output.append(outputstr)
            annots_fh.write(outputstr)
    if annots_fh is not None:
        annots_fh.close()
    outputstr = "total # of groups with shared samples:\t%d\n" % (shared_group_count)
    output.append(outputstr)
    sys.stderr.write(outputstr)
    outputstr = "total # of groups with fully annotated splices:\t%d\n" % (total_fully_annotated_count)
    output.append(outputstr)
    sys.stderr.write(outputstr)
    return output

def samples_changed(args,cache_file):
    response = urllib2.urlopen("%s/%s/samples?check_for_update=1" % (clsnapconf.SERVICE_URL,args.datasrc))
    remote_timestamp = response.read()
    remote_timestamp.rstrip()
    remote_timestamp = float(remote_timestamp)
    stats = os.stat(cache_file)
    local_timestamp = stats.st_mtime
    if remote_timestamp > local_timestamp:
        return True
    return False

def download_sample_metadata(args):
    '''Dump from Snaptron WSI the full sample metadata for a specific data compilation (source)
    to a local file if not already cached'''

    sample_records = {}
    cache_file = os.path.join(args.tmpdir,"snaptron_sample_metadata_cache.%s.tsv.gz" % args.datasrc)
    gfout = None
    if clsnapconf.CACHE_SAMPLE_METADTA:
        if os.path.exists(cache_file) and not samples_changed(args,cache_file):
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

def breakup_junction_id_query(jids):
    ln = len(jids)
    queries = []
    if ln > clsnapconf.ID_LIMIT:
        jids = list(jids)
        for i in xrange(0, ln, clsnapconf.ID_LIMIT):
            idq = 'ids='+','.join([str(z) for z in jids[i:i+clsnapconf.ID_LIMIT]])
            queries.append(idq)
    else:
        queries.append('ids='+','.join([str(z) for z in jids]))
    return queries


iterator_map = {True:SnaptronIteratorLocal, False:SnaptronIteratorHTTP}
either_patt = re.compile(r'either=(\d)')
def process_queries(args, query_params_per_region, groups, endpoint, count_function=None, local=False):
    '''General function to process the high level queries (functions) via
    one or more basic queries while tracking the results across basic queries'''

    results = {'samples':{},'queries':[],'exons':{'start':{},'end':{}},'either':0}
    if FUNCTION_TO_TYPE[args.function] == 'shared':
        results['groups_seen']={}
        results['shared']={}
        results['annotated']={}
        results['annotations']={}
    #used in intersected queries
    if FUNCTION_TO_TYPE[args.function] == 'track-queries':
        results['jids']=set()
    first = True
    group_fhs = {}
    for (group_idx, query_param_string) in enumerate(query_params_per_region):
        jids = set()
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
            if count_function is not None:
                count_function(args, results, record, group, out_fh=group_fh)
            elif args.function == INTERSECTION_FUNC:
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
        if args.function == INTERSECTION_FUNC:
            if first:
                results['jids'] = jids
            else:
                results['jids'] = results['jids'].intersection(jids)
        first = False
    #check if we have jids to query for (currently only in case of intersection function)
    if 'jids' in results and len(results['jids']) > 0:
        query_param_strings = breakup_junction_id_query(results['jids'])
        for query_param_string in query_param_strings:
            sIT = iterator_map[local](query_param_string, args.datasrc, endpoint)
            for record in sIT:
                #ignore limits and group label on this output
                sys.stdout.write(record + "\n")
    for (group,group_fh) in group_fhs.iteritems():
        group_fh.close()
    return results

    

compute_functions={PSI_FUNC:(count_samples_per_group,percent_spliced_in),JIR_FUNC:(count_samples_per_group,junction_inclusion_ratio),TRACK_EXONS_FUNC:(track_exons,filter_exons),TISSUE_SPECIFICITY_FUNC:(count_samples_per_group,tissue_specificity),SHARED_SAMPLE_COUNT_FUNC:(count_samples_per_group,report_shared_sample_counts),None:(None,None),INTERSECTION_FUNC:(None,None)}
def main(args):
    sample_records = {}
    if args.function is not None and args.function != TRACK_EXONS_FUNC:
        sample_records = download_sample_metadata(args)
    #special handling for bulk, should attempt to refactor this in the future to
    #be more streamlined
    if args.bulk_query_file is not None:
        process_bulk_queries(args)
        return
    #parse original set of queries
    (query_params_per_region, groups, endpoint) = parse_query_params(args)
    #get original functions (if passed in)
    (count_function, summary_function) = compute_functions[args.function]
    #process original queries
    results = process_queries(args, query_params_per_region, groups, endpoint, count_function=count_function, local=args.local)
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

    parser.add_argument('--function', metavar='jir', type=str, default=None, help='function to compute between specified groups of junctions ranked across samples: "jir", "psi", "ts", "ssc", and "exon"')

    parser.add_argument('--tmpdir', metavar='/path/to/tmpdir', type=str, default=clsnapconf.TMPDIR, help='path to temporary storage for downloading and manipulating junction and sample records')
    
    parser.add_argument('--psi-min-reads', metavar='PSI=mean(inclusion_junction1,inclusion_junction2)/(mean(inclusion_junction1,inclusion_junction2) + exclusion_junction)', type=int, default=clsnapconf.PSI_MIN_READS, help='minimum # of reads required in total (denominator) when calculating the PSI, default is ' + str(clsnapconf.PSI_MIN_READS))

    parser.add_argument('--limit', metavar='1', type=int, default=-1, help='# of records to return, defaults to all (-1)')

    parser.add_argument('--exon-length', metavar='50-60', type=str, default=None, help='length range of exons to look for within queried region when function is set to "exon", defaults to None (print out all exons in queried region when function is "exon")')
    parser.add_argument('--local', action='store_const', const=True, default=False, help='if running Snaptron modeules locally (skipping WSI)')
    parser.add_argument('--noraw', action='store_const', const=True, default=False, help='don\'t store output of WSI in tmpdir')
    parser.add_argument('--noheader', action='store_const', const=True, default=disable_header, help='turn off printing header in output')

    parser.add_argument('--datasrc', metavar='data_source_namd', type=str, default=clsnapconf.DS_SRAV2, help='data source instance of Snaptron to use, check clsnapconf.py for list')
    return parser


if __name__ == '__main__':
   
    parser = create_parser()
    #returned format (UCSC, and/or subselection of fields) option?
    #intersection or union of intervals option?

    args = parser.parse_args()
    if args.region is None and args.filters is None and args.metadata is None and args.query_file is None and args.bulk_query_file is None:
        sys.stderr.write("Error: no discernible arguments passed in, exiting\n")
        parser.print_help()
        sys.exit(-1)
    if args.function == TISSUE_SPECIFICITY_FUNC and 'gtex' not in args.datasrc:
        sys.stderr.write("Error: attempting to do tissue specificity (ts) on a non-GTEx Snaptron instance. Please change the SERVICE_URL setting in clsnapconf.py file to be the GTEx Snaptron instance before running this function; exiting\n")
        sys.exit(-1)
    if not os.path.exists(args.tmpdir):
        os.mkdir(args.tmpdir)
    main(args)
