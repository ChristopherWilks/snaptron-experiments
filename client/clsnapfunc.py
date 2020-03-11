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
import re
from future.utils import viewitems

import clsnapconf
import clsnaputil
import clsnapconf
from SnaptronIteratorHTTP import SnaptronIteratorHTTP
from SnaptronIteratorLocal import SnaptronIteratorLocal

GTEX_TISSUE_COL=65
#TODO: switch to using a set of classes inheriting from a base SnapFunction class
TISSUE_SPECIFICITY_FUNC='ts'
SHARED_SAMPLE_COUNT_FUNC='shared'
JIR_FUNC='jir'
TRACK_EXONS_FUNC='exon'
INTERSECTION_FUNC='intersection'
MATES_FUNC='mates'
GROUP_COV_FUNC='coverage'

apsi_region_patt = re.compile(r'^(chr[^:]+):(\d+)-(\d+)$')
def process_apsi(args, query_params_per_region, query_labels, endpoint, datasrcs, iterator_map, local=False):
    #first, don't support more than one datasrc
    if len(datasrcs) > 1:
        sys.stderr.write("APSI only supports queries on a single datasource at a time, quitting\n")
        sys.exit(-1)
    #each query params entry is really a list of query params, where the first
    rids2groups = {}
    with open(args.sample_group_file,"r") as fin:
        for line in fin:
            fields = line.rstrip().split('\t')
            group = fields[0]
            rids = fields[3].split(',')
            rids2groups.update({rid:group for rid in rids})
    groups_sorted = sorted(set(rids2groups.values()))
    sys.stdout.write("query_label\tquery\t%s\n" % ('\t'.join(groups_sorted)))
    for (qi,query_params) in enumerate(query_params_per_region):
        query_label = query_labels[qi]
        for q in query_params:
            (k,v) = q.split('=')
            params_ = []
            supporting_jx_regions = []
            if k == 'regions':
                #multiple coordinate ranges are expected
                #and assume they're in coordinate sorted order
                supporting_jx_regions = v.split(',')
                #each range is an inclusion jx, the start/end coordinate is the shared exclusion SS
                #and the end/start is the not-shared inclusion SS
            else:
                params_.append(q)
            #now create the actual queries we want to run
            params = '&'.join(params_)
            psis = []
            for (i,region) in enumerate(supporting_jx_regions):
                psis.append({})
                #init the psis per query, per group, and per numerator/denominator (n/d)
                for g in groups_sorted:
                    psis[i][g] = {'n':0.0, 'd':0.0}
                m = apsi_region_patt.search(region)
                if m is None:
                    sys.stderr.write("WARNING: region not in expected format for APSI function: %s, skipping\n" % region)
                    continue
                #main idea here: want to match a range of junctions that either:
                #1) start with the given start coordinate and have any end OR (ODD)
                #2) end with the given end coordinate and have any start (EVEN)
                chrm = m.group(1)
                start = m.group(2)
                end = m.group(3)
                either = 1
                query_param_string = "regions=%s:%s-%s&either=%d" % (chrm, start, start, either)
                inc_coord = end
                inc_col = clsnapconf.INTERVAL_END_COL
                #again, assume coordinate sorted input order
                if i % 2 != 0:
                    either = 2
                    query_param_string = "regions=%s:%s-%s&either=%d" % (chrm, end, end, either)
                    inc_coord = start
                    inc_col = clsnapconf.INTERVAL_START_COL
                if len(params_) > 0:
                    query_param_string += "&%s" % (params_)
                #run query
                #process results looking for:
                #1) the one inclusion junction which has a coordinate which matches the "inc_coord"
                #these counts will go in both the denominator AND the numerator
                #2) everything else for which the counts will go in the denominator only
                #ASSUMES: 1) 1 datasrc 2) 1 endpoint (junctions)
                #for (group_idx, query_param_string) in enumerate(query_params_per_region):
                sIT = iterator_map[local](query_param_string, datasrcs[0], endpoint)
                for record in sIT:
                    fields = record.split('\t')
                    #process record according to 
                    for sample in fields[clsnapconf.SAMPLE_IDS_COL].split(',')[1:]:
                        (rid,cov) = sample.split(':')
                        cov = int(cov)
                        if rid in rids2groups:
                            g = rids2groups[rid]
                            psis[i][g]['d'] += cov
                            if fields[inc_col] == inc_coord:
                                psis[i][g]['n'] += cov
            #out region loop but still within same event
            #now that we've processed the basic queries, we need to form the final PSIs
            #assumes only 2 basic queries for now (inclusion left/right)
            if len(psis[0]) == 0 or len(psis[1]) == 0:
                sys.stderr.write("Bad APSI query row: %s, skipping\n" % (q))
                continue
            sys.stdout.write("%s\t%s" % (query_label,q))
            for g in groups_sorted:
                #upstream psi
                psi1 = 0.0
                if g in psis[0] and psis[0][g]['d'] > 0:
                    psi1 = psis[0][g]['n'] / float(psis[0][g]['d'])
                #downstream psi
                psi2 = 0.0
                if g in psis[1] and psis[1][g]['d'] > 0:
                    psi2 = psis[1][g]['n'] / float(psis[1][g]['d'])
                psi = (psi1 + psi2)/2.0
                sys.stdout.write("\t%d" % int(round(100*psi,0)))
            sys.stdout.write('\n')


def intersect_junctions(args, results, record, group, out_fh=None):
    '''Intersect current snaptron_ids (junction ids) with previous queries results.
       If this is the last query in the group, do intersection but then print out results.
    '''
    fields = record.split('\t')
    snid = fields[clsnapconf.INTRON_ID_COL]
    return


def extract_group_coverage(args, sample_stat, group_list):
    '''Extracts the raw coverage between N junctions and adds a total sum'''
   
    covs = [0]+[sample_stat[g] for g in group_list]
    total = reduce(lambda x, y: x+y, covs)
    if total < args.min_count:
        total = -total
    covs[0] = total
    return covs

def group_coverage(args, results, group_list, sample_records):
    '''Outputs raw coverage sums for total and for each input jx group, per sample
        plus all sample metadata, one sample per line'''

    #we allow for N groups, BUT the user has
    #to specify unique group names for each individual junction
    #ahead of time e.g. L1_INCLUSION L2_INCLUSION L1_EXCLUSION, etc...
    sample_stats = results['samples']
    sample_scores = {}
    for sample in sample_stats:
        for group in group_list:
            if group not in sample_stats[sample]:
                sample_stats[sample][group]=0
        sample_scores[sample] = extract_group_coverage(args, sample_stats[sample], group_list)
    missing_sample_ids = set()
    counter = 0
    output = []
    if not args.noheader:
        rawstr = ""
        for group in group_list:
            rawstr += "%s_count\t" % group
        sheader = sample_records["header"]
        outstr = "total_count\t%s%s\n" % (rawstr,sheader)
        output.append(outstr)
        sys.stdout.write(outstr)
    for sample in sorted(sample_scores.keys(),key=lambda x: sample_scores[x][0],reverse=True):
        counter += 1
        if args.limit > -1 and counter > args.limit:
            break
        scores = sample_scores[sample]
        if sample not in sample_records:
            missing_sample_ids.add(sample)
            continue
        rawstr = "\t".join([str(x) for x in sample_scores[sample]])
        sample_record = sample_records[sample]
        outstr = "%s\t%s\n" % (rawstr,sample_record)
        output.append(outstr)
        sys.stdout.write(outstr)
    return output

def calc_psi(args, sample_stat, group_list):
    '''Calculates the simple PSI between 2 junction 
    variants where one variant has 2 junctions (inclusion),
    and the other has only one junction (exclusion).'''

    ss = sample_stat
    gl = group_list
    #TODO blatant hardcoding of group indexes here, need to allow for less rigidity (more groups)
    (inclusion1,inclusion2,exclusion) = (sample_stat[gl[0]],sample_stat[gl[1]],sample_stat[gl[2]])
    mean_inclusion = (inclusion1 + inclusion2) / 2.0
    total = mean_inclusion + exclusion
    psi = mean_inclusion / float(total)
    if inclusion1 == 0 or inclusion2 == 0 or total < args.min_count:
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
    majority_psi_sample_count = 0
    exclusion_sample_count = 0
    for sample in sample_stats:
        for group in group_list:
            if group not in sample_stats[sample]:
                sample_stats[sample][group]=0
        sample_scores[sample] = calc_psi(args, sample_stats[sample], group_list)
        if args.summarize:
            (psi,total,inclusion1,inclusion2,exclusion) = sample_scores[sample]
            if psi > 0.5:
                majority_psi_sample_count+=1
            if exclusion > 0:
                exclusion_sample_count+=1
    missing_sample_ids = set()
    counter = 0
    output = []
    if args.summarize:
        outstr = '%s: # samples with > 0.5 PSI: %d, #samples with an exclusion jx %d, ratio: %s\n' % (args.datasrc,majority_psi_sample_count, exclusion_sample_count, majority_psi_sample_count/float(exclusion_sample_count))
        sys.stderr.write(outstr)
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


def junction_inclusion_ratio(args, results, group_list, sample_records, print_output=True):
    '''Calculates the junction inclusion ratio across basic queries for
    for the groups of junctions specified in group_list for samples
    from sample_records; junction group coverages are reported in results per sample'''

    no_sample_records = False
    if sample_records is None:
        sample_records = {}
        sample_records["header"] = ""
        no_sample_records = True

    sample_stats = results['samples']
    group_a = group_list[0]
    group_b = group_list[1]
    sample_scores = {}
    A_majority_jir_sample_count = 0
    B_majority_jir_sample_count = 0
    neither_sample_count = 0
    for sample in sample_stats:
        if group_a not in sample_stats[sample]:
            sample_stats[sample][group_a]=0
        if group_b not in sample_stats[sample]:
            sample_stats[sample][group_b]=0
        #need to have a minimum number of reads in the denominator
        if (sample_stats[sample][group_a] + sample_stats[sample][group_b]) < args.min_count_jir:
            sample_scores[sample] = -100
        else:
            sample_scores[sample] = calc_jir(sample_stats[sample][group_a], sample_stats[sample][group_b])
        if args.summarize:
            if sample_scores[sample] > 0:
                B_majority_jir_sample_count+=1
            elif sample_scores[sample] < 0:
                A_majority_jir_sample_count+=1
            else:
                neither_sample_count+=1
    if args.summarize:
        total_sample_count = A_majority_jir_sample_count + B_majority_jir_sample_count + neither_sample_count
        outstr = '%s: # samples > 0 JIR: %d, # samples < 0 jir: %d, # samples == 0: %d, ratios: >0/total: %s, <0/total: %s\n' % (args.datasrc, B_majority_jir_sample_count, A_majority_jir_sample_count, neither_sample_count, B_majority_jir_sample_count/float(total_sample_count), A_majority_jir_sample_count/float(total_sample_count))
        sys.stderr.write(outstr)
    missing_sample_ids = set()
    counter = 0
    output = []
    if not args.noheader:
        sheader = sample_records["header"]
        outstr = "jir_score\t%s raw count\t%s raw count\t%s\n" % (group_a,group_b,sheader)
        output.append(outstr)
        sys.stdout.write(outstr)
    for sample in sorted(sample_scores.keys(),key=sample_scores.__getitem__,reverse=True):
        score = sample_scores[sample]
        #should be at the bottom of the samples where everything else is disqualified
        if score == -100:
            break
        counter += 1
        if args.limit > -1 and counter > args.limit:
            break
        if no_sample_records:
            sample_records[sample]=""
        elif sample not in sample_records:
            missing_sample_ids.add(sample)
            continue
        sample_record = sample_records[sample]
        outstr = "%s\t%d\t%d\t%s\n" % (str(score),sample_stats[sample][group_a],sample_stats[sample][group_b],sample_record)
        output.append(outstr)
        sys.stdout.write(outstr)
    return output

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


FUNCTION_TO_TYPE={MATES_FUNC:'not-shared', TRACK_EXONS_FUNC:'not-shared', JIR_FUNC:'not-shared', clsnaputil.PSI_FUNC:'not-shared', None:None, TISSUE_SPECIFICITY_FUNC:'shared',SHARED_SAMPLE_COUNT_FUNC:'shared',INTERSECTION_FUNC:'track-queries',GROUP_COV_FUNC:'not-shared'}
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

region_patt = re.compile(r'^(chr[0-9a-zA-Z_\-]+):(\d+)-(\d+)$')
def report_splice_mates(args, results, group_list, sample_records):
    '''Takes donor (acceptor) and finds all junctions originating
        from that site.  It also finds the coverage of the base
        just beyond it in the intron in all samples.  It then
        computes the splice event score for every sample for
        each junction. This is the splice mate "reporter" function 
        for the more general sum_sample_coverage "compute" function.'''
    subdelim='\t'
    totals = results['all_sample_sums']
    #do RI point query
    #curl "http://snaptron.cs.jhu.edu/supermouse/bases?regions=chr4:156177822-156178635"
    #Null out params which won't be used by the call to the bases endpoint
    args.filters=None
    args.metadata=None
    #force match to inclusive base 1 end
    m = region_patt.search(args.region)
    chrom = m.group(1)
    start = int(m.group(2))
    end = int(m.group(3))
    args.region = '%s:%s-%s' % (chrom,str(start),str(end))
    (query_param_strings, groups, endpoint, datasrcs) = clsnaputil.parse_command_line_args(args)
    use_local = False
    #assume one query
    #for now, only support a single datasrc via this interface
    if len(datasrcs) != 1:
        sys.stderr.write("# of datasrc/compilations != 1: %s, this function only supports one datasrc/compilation, quiting\n" % (','.join(datasrcs)))
        sys.exit(-1)
    sIT = results['siterator'](query_param_strings[0], args.datasrc, clsnapconf.BASES_ENDPOINT)
    #expecting only 2 rows for this point query, the header and the base itself
    sample_ids = sIT.next().split('\t')[clsnapconf.INTERVAL_END_COL:]
    base_vals = sIT.next().split('\t')
    base_range_fields = base_vals[:clsnapconf.INTERVAL_END_COL]
    #a little hack to ensure the bases results line up with the junction columns
    base_range_fields[0]="%s\t-1" % base_range_fields[0]
    base_vals = [float(base_val) for base_val in base_vals[clsnapconf.INTERVAL_END_COL:]]
    for i, base_val in enumerate(base_vals):
        if base_val == 0:
            continue
        if sample_ids[i] not in totals:
            totals[sample_ids[i]] = 0
        totals[sample_ids[i]] += base_val
    
    samples = sorted(totals.keys(), key=int)
    sample_header = subdelim.join(samples)
    sys.stdout.write("\t".join(results['header_fields']) + "\t" + sample_header+"\n")

    #filter out low values from denominator (totals), but keep for final output as 0's
    totals = {sid:float(x) for sid,x in viewitems(totals) if float(x) >= args.min_count}

    #print out RI counts
    sid_off = clsnapconf.BASE_SAMPLE_ID_OFFSETS[args.datasrc]
    all_samples = {x:0 for x in samples}
    all_samples.update({str(sample_ids[i]):(x/totals[str(sample_ids[i])]) for i,x in enumerate(base_vals) if x != 0 and str(sample_ids[i]) in totals})
    filler_length = (clsnapconf.SAMPLE_IDS_COL - clsnapconf.INTERVAL_END_COL)
    sys.stdout.write("\t".join(base_range_fields)+('\t'*filler_length))
    sys.stdout.write(subdelim.join([str(all_samples[s]) for s in samples])+"\n")

    #now calculate mate score for all junctions if an event type is not specified
    if not args.event_type:
        junctions = results['junctions']
        for jx_id in sorted(junctions.keys(), key=int):
            jx = junctions[jx_id]
            all_samples = {x:0 for x in samples}
            all_samples.update({x.split(":")[0]:(int(x.split(":")[1])/float(totals[x.split(":")[0]])) for x in jx[clsnapconf.SAMPLE_IDS_COL].split(',')[1:] if x.split(":")[0] in totals})
            sys.stdout.write("\t".join(jx[:clsnapconf.SAMPLE_IDS_COL])+"\t")
            sys.stdout.write(subdelim.join([str(all_samples[s]) for s in samples])+"\n")


def sum_sample_coverage(args, results, record, group, out_fh=None):
    '''Sums coverage for every splice junction for every sample
        across all junctions while also tracking per-sample coverage
        for each junction separately.'''

    if 'all_sample_sums' not in results:
        results['all_sample_sums'] = {}
    if 'junctions' not in results:
        results['junctions'] = {}
    if out_fh is not None:
        out_fh.write(record+"\n")
    fields = record.split('\t')
    if 'snaptron_id' == fields[clsnapconf.INTRON_ID_COL]:
        results['header_fields']=fields[:clsnapconf.SAMPLE_IDS_COL]
        return
    samples = fields[clsnapconf.SAMPLE_IDS_COL][1:].split(',')
    #results['sample_covs'][fields[clsnapconf.INTRON_ID_COL]] = [x.split(":") for x in samples]
    results['junctions'][fields[clsnapconf.INTRON_ID_COL]] = fields
    sample_covs = [x.split(":")[1] for x in samples]
    samples = [x.split(":")[0] for x in samples]
    #need to track counts for the denominator
    for i,sid in enumerate(samples):
        if sid not in results['all_sample_sums']:
            results['all_sample_sums'][sid]=0
        results['all_sample_sums'][sid]+=int(sample_covs[i])


