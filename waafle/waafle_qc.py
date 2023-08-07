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

from waafle import utils as wu

# ---------------------------------------------------------------
# description
# ---------------------------------------------------------------

description = wu.describe( """
{SCRIPT}: Applies junction results to QC WAAFLE calls.

Filter WAAFLE lgt output to require that junctions
be supported by mate-pair or coverage evidence.
""" )

# ---------------------------------------------------------------
# output formats
# ---------------------------------------------------------------

c_formats = {}
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
        "contig_profile",
        help="lgt output from waafle_orgscorer (tsv format)",
        )
    g.add_argument( 
        "junctions",
        help="output from waafle_junctions for contigs of interest",
        )

    g = parser.add_argument_group( "filtering parameters" )
    g.add_argument( 
        "--min-junction-hits",
        type=int,
        default=2,
        metavar="<int>",
        help="minimum read-hits to 'ok' a junction\n[default: 2]",
        )
    g.add_argument( 
        "--min-junction-ratio",
        type=int,
        default=0.5,
        metavar="<float>",
        help="minimum coverage (relative to flanking genes) to 'ok' a junction\n[default: 0.5]",
        )
    """
    g.add_argument( 
        "--min-contig-genes",
        type=int,
        default=2,
        metavar="<int>",
        help="minimum gene count for the contig\n[default: 1]",
        )
    g.add_argument( 
        "--min-contig-length",
        type=int,
        default=500,
        metavar="<int>",
        help="minimum length for the contig\n[default: 500]",
        )
    """

    g = parser.add_argument_group( "misc options" )
    """
    g.add_argument( 
        "--ab-only",
        action="store_true",
        help="Only apply filtering to LGT junctions (AB/BA)",
        )
    """
    g.add_argument( 
        "--outfile",
        type=str,
        default=None,
        metavar="<path>",
        help="Path for filtered outputs\n[default: derive from input]",
        )

    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):

    args = get_args( )

    # load junctions data
    hits = {}
    covs = {}
    wu.say( "Loading junctions report." )
    F = wu.Frame( args.junctions )
    # loop over junctions
    for R in F.iter_rowdicts( ):
        contig = R["CONTIG"]
        gene1 = R["GENE1"]
        gene2 = R["GENE2"]
        hits.setdefault( contig, {} )[(gene1, gene2)] = int( R["JUNCTION_HITS"] )
        covs.setdefault( contig, {} )[(gene1, gene2)] = float( R["RATIO"] )

    # filter contigs
    total = 0
    failed = 0
    outfile = args.outfile
    if outfile is None:
        outfile = args.contig_profile + ".qc_pass"
    # load results, open new file, write headers
    F = wu.Frame( args.contig_profile )
    fh = wu.try_open( outfile, "w" )
    wu.write_rowdict( None, F.headers, file=fh )
    # loop over contigs
    for R in F.iter_rowdicts( ):
        total += 1
        contig = R["CONTIG_NAME"]
        # contig-level filters
        if contig not in hits or contig not in covs:
            failed += 1
            wu.say( "Missing junction data for contig:", contig )
            continue        
        loci = R["LOCI"].split( "|" )
        synteny = R["SYNTENY"]
        qc_pass = True
        for i in range( len( loci ) - 1 ):
            my_test = True
            spair = synteny[i] + synteny[i+1]
            if spair not in ["AB", "BA"]:
                continue
            gpair = (loci[i], loci[i+1])
            my_hits = hits[contig].get( gpair, -1 )
            my_hits = my_hits >= args.min_junction_hits
            my_covs = covs[contig].get( gpair, -1 )
            my_covs = my_covs >= args.min_junction_ratio
            my_test = my_hits or my_covs
            qc_pass = qc_pass and my_test
        if not qc_pass:
            failed += 1
            wu.say( "Failed QC:", contig )
        else:
            wu.write_rowdict( R, F.headers, file=fh )

    # wrap-up
    wu.say( "Failure rate: {} of {} ({:.1f}%)".format( failed, total, 100 * failed / float( total ) ) )
    wu.say( "Finished successfully." )           
                
if __name__ == "__main__":
    main( )
