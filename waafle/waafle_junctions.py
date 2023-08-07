#!/usr/bin/env python

"""
This module is a part of:
WAAFLE, a [W]orkflow to [A]nnotate [A]ssemblies and [F]ind [L]GT [E]vents

Copyright (c) 2023 Harvard T.H. Chan School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from __future__ import print_function # Python 2.7+ required
import os
import sys
import csv
import argparse
import re
from collections import Counter

import numpy as np

from waafle import utils as wu

# ---------------------------------------------------------------
# description
# ---------------------------------------------------------------

description = wu.describe( """
{SCRIPT}: Generate gene-gene junction stats for contig QC

This script maps reads to contigs (or analyzes an existing mapping)
to find support for gene-gene junctions. A junction is supported if
there are mate-pairs that span the junction. As this may be infeasible
for larger junctions (e.g. those exceeding a typical insert size of 300 nt),
a junction can also be supported by having a reasonable average coverage
relative to its flanking genes.
""" )

# ---------------------------------------------------------------
# output formats
# ---------------------------------------------------------------

c_formats = {}

c_formats["site_hits"] = """
contig
mean
stdev
depths
"""

c_formats["gene_hits"] = """
contig
gene1
gene2
hits
"""

c_formats["junctions"] = """
contig
gene1
gene2
len_gene1
len_gene2
gap
junction_hits
coverage_gene1
coverage_gene2
coverage_junction
ratio
"""

for name, items in c_formats.items( ):
    c_formats[name] = [k for k in items.split( "\n" ) if k != ""]

# ---------------------------------------------------------------
# cli
# ---------------------------------------------------------------

def get_args( ):

    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
        )

    g = parser.add_argument_group( "required inputs" )
    g.add_argument( 
        "contigs",
        help="contigs file (fasta format)",
        )
    g.add_argument( 
        "gff",
        help="GFF file for provided contigs",
        )

    g = parser.add_argument_group( "provide paired reads or a .sam file" )
    g.add_argument( 
        "--reads1",
        metavar="<path>",
        help="sequencing reads (mate-1)",
        )
    g.add_argument( 
        "--reads2",
        metavar="<path>",
        help="sequencing reads (mate-2)",
        )
    g.add_argument( 
        "--sam",
        metavar="<path>",
        help="sam file (from existing alignment)",
        )

    g = parser.add_argument_group( "output options" )
    g.add_argument( 
        "--tmpdir",
        default=".",
        metavar="<path>",
        help="where to place temp outputs\n[default: ./]",
        )
    g.add_argument( 
        "--outdir",
        default=".",
        metavar="<path>",
        help="where to place main outputs\n[default: ./]",
        )
    g.add_argument( 
        "--basename",
        metavar="<str>",
        help="basename for output files\n[default: <derived from input>]",
        )
    g.add_argument( 
        "--write-detailed-output",
        action="store_true",
        help="write out coverage values for all sites and all junctions\n[default: off]",
        )

    g = parser.add_argument_group( "filtering parameters" )
    g.add_argument(
        "--min-overlap-sites",
        type=int,
        default=25,
        metavar="<int>",
        help="minimum nucleotide overlap for counting a read-gene hit\n[default: 25]",
        )

    g = parser.add_argument_group( "bowtie2 options" )
    g.add_argument( 
        "--bowtie2-build",
        default="bowtie2-build",
        metavar="<path>",
        help="path to bowtie2-build\n[default: $PATH]",
        )
    g.add_argument( 
        "--bowtie2",
        default="bowtie2",
        metavar="<path>",
        help="path to bowtie2\n[default: $PATH]",
        )
    g.add_argument( 
        "--threads",
        type=int,
        default=1,
        metavar="<int>",
        help="number of threads for bowtie2 steps\n[default: 1]",
        )
    g.add_argument( 
        "--resume",
        action="store_true",
        help="if set, use existing .index and/or .sam if found\n[default: off]",
        )

    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utils for running bowtie2
# ---------------------------------------------------------------

def bowtie2_build( p_bowtie2_build=None, p_contigs=None, p_index=None, 
                   args=None, ):
    alias = {
        "PROG":    p_bowtie2_build,
        "CONTIGS": p_contigs, 
        "INDEX":   p_index, 
        }
    if args.resume and os.path.exists( p_index + ".1.bt2" ):
        wu.say( "RESUMING: The index <{INDEX}> already exists.".format( **alias ) )
    else:
        wu.say( "Indexing <{CONTIGS}> to <{INDEX}>.".format( **alias ) )
        command = [
            "{PROG}",
            "{CONTIGS}",
            "{INDEX}",
            ]
        command = " ".join( command )
        command = command.format( **alias )
        os.system( command )
        wu.say( "Build complete." )
    return None

def bowtie2_align( p_bowtie2=None, p_reads1=None, p_reads2=None, 
                   p_index=None, p_sam=None, args=None, ):
    alias = {
        "PROG":    p_bowtie2,
        "READS1":  p_reads1,
        "READS2":  p_reads2,
        "INDEX":   p_index,
        "SAM":     p_sam,
        "THREADS": args.threads,
        }
    if args.resume and os.path.exists( p_sam ):
        wu.say( "RESUMING: A sam mapping <{SAM}> already exists.".format( **alias ) )
    else:
        wu.say( "Performing bowtie2 alignment." )
        command = [
            "{PROG}",
            "-x {INDEX}",
            "-1 {READS1}",
            "-2 {READS2}",
            "-S {SAM}",
            "--threads {THREADS}",
            "--no-mixed",
            "--no-discordant",
            ]
        command = " ".join( command )
        command = command.format( **alias )
        os.system( command )
        wu.say( "Alignment complete." )
    return None

# ---------------------------------------------------------------
# utils for parsing SAM/GFF comparison
# ---------------------------------------------------------------

def concordant_hits( p_sam=None, ):
    counter = 0
    mate1 = None
    mate2 = None
    for hit in wu.iter_sam_hits( p_sam ):
        # progress
        counter += 1
        if counter % int( 1e5 ) == 0:
            wu.say( "  SAM alignments processed: {:.1f}M".format( counter / 1e6 ) )
        # weave
        mate1 = mate2
        mate2 = hit
        # edge case
        if mate1 is None:
            continue
        # not a mate pair
        elif mate1.qseqid != mate2.qseqid:
            continue
        # not concordant
        elif mate1.sseqid != mate2.sseqid:
            continue
        # good pair
        else:
            yield [mate1, mate2]

def find_hit_loci( mate1=None, mate2=None, loci=None, args=None ):
    hits = set( )
    for L in loci:
        a1, a2 = L.start, L.end
        for read in [mate1, mate2]:
            b1, b2 = read.sstart, read.send
            overlap = wu.calc_overlap( a1, a2, b1, b2, normalize=False )
            if overlap >= args.min_overlap_sites:
                hits.add( L.code )
    return hits

# ---------------------------------------------------------------
# utils for evaluating a contig
# ---------------------------------------------------------------

def evaluate_contig( loci=None, coverage=None, gene_hits=None, args=None, ):
    rowdicts = []
    loci = sorted( loci, key=lambda x: x.start )
    for i in range( len( loci ) - 1 ):
        rowdict = {}
        L1 = loci[i]
        L2 = loci[i+1]
        rowdict["gene1"]             = my_code1 = L1.code
        rowdict["gene2"]             = my_code2 = L2.code
        rowdict["len_gene1"]         = len( L1 )
        rowdict["len_gene2"]         = len( L2 )
        rowdict["gap"]               = my_gap = L2.start - L1.end - 1
        # check hits
        rowdict["junction_hits"]     = my_hits = gene_hits.get( (my_code1, my_code2), 0 )
        # check coverage (note: base-0 start and pythonic end)
        rowdict["coverage_gene1"]    = my_cov1 = np.mean( coverage[L1.start-1:L1.end] )
        rowdict["coverage_gene2"]    = my_cov2 = np.mean( coverage[L2.start-1:L2.end] )
        # define junction coverage as 0 if there's no gap
        my_covj = 0.0 if my_gap <= 0 else np.mean( coverage[L1.end-1:L2.start] )
        rowdict["coverage_junction"] = my_covj
        # note: pseudocount to avoid /0
        my_mean = np.mean( [my_cov1, my_cov2] )
        rowdict["ratio"]             = my_ratio = my_covj / ( my_mean + 1e-6 )
        rowdicts.append( rowdict )
    return rowdicts

# ---------------------------------------------------------------
# utils for output formatting
# ---------------------------------------------------------------

def write_detailed_output( basename=None, outdir=None,
                           contig_coverage=None, contig_hits=None, ):

    # first file can get pretty big, hence gzip
    p_site_hits = wu.name2path( basename, outdir, ".site_hits.tsv.gz" )
    p_gene_hits = wu.name2path( basename, outdir, ".gene_hits.tsv" )

    # write: site_hits
    wu.say( "Writing site hits." )
    with wu.try_open( p_site_hits, "w" ) as fh:
        wu.write_rowdict( 
            format=c_formats["site_hits"], 
            file=fh, )
        for c in sorted( contig_coverage ):
            depths = contig_coverage[c]
            rowdict = {
                "contig" : c,
                "mean"   : np.mean( depths ),
                "stdev"  : np.std( depths ),
                "depths" : " ".join( ["{:.0f}".format( k ) for k in depths] ),
                }
            wu.write_rowdict( 
                rowdict=rowdict, 
                format=c_formats["site_hits"], 
                file=fh, )

    # write: gene_hits
    wu.say( "Writing gene-pair hits." )
    with wu.try_open( p_gene_hits, "w" ) as fh:
        wu.write_rowdict( 
            format=c_formats["gene_hits"], 
            file=fh, )
        for c in sorted( contig_hits ):
            for code1, code2 in sorted( contig_hits[c] ):
                if code2 > code1:
                    continue
                value = contig_hits[c][(code1, code2)]
                rowdict = {
                    "contig" : c,
                    "gene1"  : code1,
                    "gene2"  : code2,
                    "hits"   : value,
                    }
                wu.write_rowdict( 
                    rowdict=rowdict, 
                    format=c_formats["gene_hits"], 
                    file=fh, )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):

    # begin
    args = get_args( )
    p_contigs = args.contigs
    p_gff = args.gff

    # define files
    p_outdir = args.outdir
    p_tmpdir = args.tmpdir
    basename = args.basename
    if basename is None:
        basename = wu.path2name( p_contigs )
    p_index     = wu.name2path( basename, p_tmpdir, ".index" )
    p_sam       = wu.name2path( basename, p_tmpdir, ".sam" )
    p_junctions = wu.name2path( basename, p_outdir, ".junctions.tsv" )

    # alignment workflow
    if args.sam is not None:
        p_sam = args.sam
        wu.say( "Using specified SAM file:", p_sam )
    elif args.reads1 is not None and args.reads2 is not None:
        # build process
        bowtie2_build( 
            p_bowtie2_build=args.bowtie2_build,
            p_contigs=args.contigs,
            p_index=p_index,
            args=args,
            )
        # alignment process
        bowtie2_align(
            p_bowtie2=args.bowtie2,
            p_reads1=args.reads1, 
            p_reads2=args.reads2,
            p_index=p_index,
            p_sam=p_sam,
            args=args,
            )
    else:
        wu.die( "Must provide READS or SAM file." )

    # load contig data
    wu.say( "Loading contig lengths." )
    contig_lengths = wu.read_contig_lengths( p_contigs )
    contig_coverage = {}
    for name, length in contig_lengths.items( ):
        contig_coverage[name] = np.zeros( length )
    wu.say( "Loading contig gene coordinates." )
    contig_loci = {}
    for name, loci in wu.iter_contig_loci( p_gff ):
        contig_loci[name] = loci
    contig_hits = {}

    # post-processing workflow
    wu.say( "Processing SAM file." )
    for mate1, mate2 in concordant_hits( p_sam ):
        contig = mate1.sseqid
        inner = contig_hits.setdefault( contig, Counter( ) )
        # update pers-site coverage (note: base-0 start and pythonic end)
        coords = [mate1.sstart, mate1.send, mate2.sstart, mate2.send]
        L = min( coords ) - 1
        R = max( coords ) - 1
        contig_coverage[contig][L:R+1] += 1
        # find hit loci
        hits = find_hit_loci( 
            mate1=mate1, 
            mate2=mate2, 
            loci=contig_loci.get( contig, [] ),
            args=args,
            )
        # attach self counts
        for code in hits:
            inner[(code, code)] += 1
        # attach pair counts (note: symmetric storage for safer lookup)
        for code1 in hits:
            for code2 in hits:
                if code1 != code2:
                    inner[(code1, code2)] += 1

    # detailed output?
    if args.write_detailed_output:
        write_detailed_output(
            basename=basename,
            outdir=p_outdir,
            contig_coverage=contig_coverage,
            contig_hits=contig_hits,
            )

    # write junction report
    wu.say( "Writing junction report." )
    with wu.try_open( p_junctions, "w" ) as fh:
        wu.write_rowdict( 
            format=c_formats["junctions"], 
            file=fh, )
        for c in sorted( contig_lengths ):
            rowdicts = evaluate_contig( 
                loci=contig_loci.get( c, [] ),
                coverage=contig_coverage[c], 
                gene_hits=contig_hits.get( c, {} ),
                args=args,
                )
            for rowdict in rowdicts:
                rowdict["contig"] = c
                wu.write_rowdict( 
                    rowdict=rowdict, 
                    format=c_formats["junctions"], 
                    file=fh, )

    # end
    wu.say( "Finished successfully." )

if __name__ == "__main__":
    main( )
