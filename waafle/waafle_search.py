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
{SCRIPT}: Step 1 in the WAAFLE pipeline

This script executes a custom BLAST search of a set of contigs against 
a WAAFLE-formatted database.
""" )

# ---------------------------------------------------------------
# cli
# ---------------------------------------------------------------

def get_args( ):
    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
        )
    parser.add_argument( 
        "query",
        help="contigs file (fasta format)",
        )
    parser.add_argument( 
        "db",
        help="path to WAAFLE BLAST database",
        )
    parser.add_argument( 
        "--blastn",
        default="blastn",
        metavar="<path>",
        help="path to blastn binary\n[default: $PATH]",
        )
    parser.add_argument( 
        "--threads",
        default="1",
        metavar="<int>",
        help="number of CPU cores to use in blastn search\n[default: 1]",
        )
    parser.add_argument( 
        "--out",
        default=None,
        metavar="<path>",
        help="path for blast output file\n[default: <derived from input>]",
        )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):
    args = get_args( )
    if args.out is None:
        name = wu.path2name( args.query )
        args.out = wu.name2path( name, ".", ".blastout" )
    params = {
        "QUERY"   : args.query,
        "OUTFILE" : args.out,
        "BLASTN"  : args.blastn,
        "DB"      : args.db,
        "MAXTAR"  : wu.c_max_target_seqs,
        "THREADS" : args.threads,
        "FORMAT"  : wu.c_blast_format_string,
        }
    command = [
        "{BLASTN}",
        "-query {QUERY}",
        "-db {DB}",
        "-out {OUTFILE}",
        "-max_target_seqs {MAXTAR}",
        "-num_threads {THREADS}",
        "-outfmt \'{FORMAT}\'",
        ]
    command = " ".join( command ).format( **params )
    wu.run_process(command)

if __name__ == "__main__":
    main( )
