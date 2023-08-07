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

from waafle import utils as wu

# ---------------------------------------------------------------
# description
# ---------------------------------------------------------------

description = wu.describe( """
{SCRIPT}: (Optional) Step 1.5 in the WAAFLE pipeline

Use the results of waafle_search to identify candidate gene
loci in a set of contigs and output them as a GFF file for use
in the next step. Users can optionally supply their own (independently-generated)
GFF file.
""" )

# ---------------------------------------------------------------
# cli
# ---------------------------------------------------------------

def get_args( ):
    """ Get arguments passed to script """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument( 
        "blastout",
        help="(custom) blast output from waafle_search",
        )
    parser.add_argument( 
        "--gff",
        default=None,
        metavar="<path>",
        help="path for (output) waafle gene calls (.gff)\n[default: <derived from input>]",
        )
    parser.add_argument(
        "--min-overlap",
        default=0.1,
        type=float,
        metavar="<float>",
        help=("if a large hit covers this fraction of a smaller hit, "
              "consider them part of the same gene group\n[default: 0.1]"),
        )
    attach_shared_args( parser )
    args = parser.parse_args( )
    return args

def attach_shared_args( parser ):
    """ these arguments are shared with orgscorer """
    parser.add_argument(
        "--min-gene-length",
        default=200,
        type=float,
        metavar="<int>",
        help="minimum allowed gene length\n[default: 200]",
        )
    parser.add_argument(
        "--min-scov",
        default=0.75,
        type=float,
        metavar="<float>",
        help="(modified) scoverage filter for hits to gene catalog\n[default: 0.75]",
        )
    parser.add_argument(
        "--stranded",
        action="store_true",
        help="only merge hits into hits/genes of the same strandedness\n[default: off]",
        )

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def hits2ints( hits, scov ):
    """ filter hits by scoverage; convert to intervals """
    intervals = []
    for hit in hits:
        if hit.scov_modified >= scov:
            intervals.append( [hit.qstart, hit.qend, hit.sstrand] )
    return intervals

def overlap_inodes( inode1, inode2 ):
    """ compute overlap between two intervals """
    a1, b1 = inode1.start, inode1.stop
    a2, b2 = inode2.start, inode2.stop
    overlap = wu.calc_overlap( a1, b1, a2, b2 )
    return overlap
    
def merge_inodes( *inodes ):
    """ merge overlapping intervals into a single node """
    start = stop = None
    strand_rank = []
    for inode in inodes:
        if start is None or inode.start < start:
            start = inode.start
        if stop is None or inode.stop > stop:
            stop = inode.stop
        strand_rank.append( [len( inode ), inode.strand] )
    # assign strand of largest interval to the merged interval
    strand = sorted( strand_rank )[-1][1]
    return wu.INode( start, stop, strand )

def make_inodes( intervals ):
    """ convert list of intervals to list of inodes """
    inodes = []
    for interval in intervals:
        start, stop = interval[0:2]
        strand = "+" if len( interval ) < 3 else interval[2]
        inodes.append( wu.INode( start, stop, strand ) )
    return inodes
        
def overlap_intervals( intervals, threshold, stranded ):
    """ find and collapse overlapping intervals """
    inodes = make_inodes( intervals )
    inodes = sorted( inodes, key=lambda inode: inode.start )
    for index, inode1 in enumerate( inodes ):
        for inode2 in inodes[index+1:]:
            if not stranded or inode1.strand == inode2.strand:
                score = overlap_inodes( inode1, inode2 )
                if score >= threshold:
                    # store as edge in network
                    inode1.attach( inode2 )
                    inode2.attach( inode1 )
                elif score == 0:
                    # no further inode2 can overlap this inode1
                    break
    # divide intervals into "connected components"
    ccs = []
    for inode in inodes:
        if not inode.visited:
            ccs.append( inode.get_connected_component( ) )
    # merge intervals and report as simple lists
    new_intervals = []
    for cc in ccs:
        merged = merge_inodes( *cc ).to_list( )
        new_intervals.append( merged )
    return new_intervals

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

"""
NOTES:

1) Organize BLAST hits that pass the scov filter into intervals
2) Sort intervals by start site and find connected (overlapping) intervals
3) Merge connected intervals into genes, and filter by gene length
4) Print genes as GFF

GFF FORMAT:

1. seqname
2. source
3. feature
4. start
5. end
6. score
7. strand
9. frame
8. attribute
"""

def main( ):

    args = get_args( )
    if args.gff is None:
        name = wu.path2name( args.blastout )
        args.gff = wu.name2path( name, ".", ".gff" )

    fh_gff = wu.try_open( args.gff, "w" )
    writer = csv.writer( fh_gff, csv.excel_tab )

    for contig, hits in wu.iter_contig_hits( args.blastout ):
        intervals = hits2ints( 
            hits, 
            args.min_scov,
            )
        intervals = overlap_intervals( 
            intervals, 
            args.min_overlap, 
            args.stranded == "on",
            )
        for start, stop, strand in intervals:
            gene_length = stop - start + 1
            if gene_length >= args.min_gene_length:
                items = [
                    contig,
                    "waafle_genecaller",
                    "gene",
                    start,
                    stop,
                    ".",
                    strand,
                    0,
                    ".",
                    ]
                writer.writerow( [str( k ) for k in items] )

    fh_gff.close( )
    wu.say( "Finished successfully." )

if __name__ == "__main__":
    main( )
