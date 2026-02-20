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

from __future__ import print_function  # Python 2.7+ required

import os
import sys
import re
import gzip
import bz2
import io
import subprocess
from collections import OrderedDict

import numpy as np

PY3 = (sys.version_info[0] >= 3)

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# GENERIC HELPER CLASSES / FUNCTIONS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

def run_process(command, success_message="\nFinished successfully."):
    issues = 0
    say("Executing command:", command)
    try:
        stdout = subprocess.check_output(command, shell=True, stderr=subprocess.STDOUT)
        if not isinstance(stdout, str):
            stdout = stdout.decode("utf-8", errors="replace")
        say(stdout)
        low = stdout.lower()
        if "warning" in low:
            say("\nFinished with warnings.")
            issues += 1
        elif "error" in low:
            say("\nFinished with errors.")
            issues += 1
        else:
            say(success_message)
    except (subprocess.CalledProcessError, EnvironmentError):
        say("\nFinished with errors.")
        issues += 1
    return issues


def say(*args):
    print(" ".join(map(str, args)), file=sys.stderr)


def die(*args):
    args = ["LETHAL ERROR:"] + list(args)
    say(*args)
    sys.exit("EXITING.")


def path2name(p):
    return os.path.split(p)[1].split(".")[0]


def name2path(name, root=".", ext=""):
    return os.path.join(root, name + ext)


# ---------------------------------------------------------------
# File open helpers
# ---------------------------------------------------------------

def _bz2_open(path, mode):
    # Python 2 compatibility: bz2.open may not exist
    if hasattr(bz2, "open"):
        return bz2.open(path, mode)
    return bz2.BZ2File(path, mode)


def try_open_bin(path, mode="rb", buffering=1024 * 1024):
    """
    Fast binary open with big buffering for plain/gz/bz2.
    Use mode "rb" for reading, "wb" for writing.
    """
    try:
        if path.endswith(".gz"):
            raw = gzip.open(path, mode) if hasattr(gzip, "open") else gzip.GzipFile(path, mode)
            if "r" in mode:
                return io.BufferedReader(raw, buffer_size=buffering)
            return raw
        elif path.endswith(".bz2"):
            raw = _bz2_open(path, mode)
            if "r" in mode:
                return io.BufferedReader(raw, buffer_size=buffering)
            return raw
        else:
            # buffering arg only supported in Py3 open() signature
            return open(path, mode, buffering=buffering) if PY3 else open(path, mode)
    except Exception as e:
        sys.exit("Can't open file: {} ({})".format(path, e))


def try_open_text(path, mode="rt", encoding="utf-8", errors="replace", newline=""):
    """
    Text open for plain/gz/bz2 safe for print(..., file=fh).
    """
    try:
        if path.endswith(".gz"):
            raw = gzip.open(path, mode.replace("t", "").replace("b", "") + "b") if hasattr(gzip, "open") else gzip.GzipFile(path, mode.replace("t", "").replace("b", "") + "b")
            if PY3:
                return io.TextIOWrapper(raw, encoding=encoding, errors=errors, newline=newline)
            return io.TextIOWrapper(raw, encoding=encoding, errors=errors)

        elif path.endswith(".bz2"):
            raw = _bz2_open(path, mode.replace("t", "").replace("b", "") + "b")
            if PY3:
                return io.TextIOWrapper(raw, encoding=encoding, errors=errors, newline=newline)
            return io.TextIOWrapper(raw, encoding=encoding, errors=errors)

        else:
            if PY3:
                return open(path, mode, encoding=encoding, errors=errors, newline=newline)
            return io.open(path, mode, encoding=encoding, errors=errors)
    except Exception as e:
        sys.exit("Can't open file: {} ({})".format(path, e))


# Legacy wrapper (kept for compatibility with other WAAFLE modules)
def try_open(path, *args):
    mode = args[0] if len(args) >= 1 else "rt"
    if "b" in mode:
        return try_open_bin(path, mode=mode)
    return try_open_text(path, mode=mode)


# ---------------------------------------------------------------
# Fast tab parsing
# ---------------------------------------------------------------

def iter_tab_rows_bytes(path, comment_prefix=None):
    """
    Yield rows as list[bytes] split on TAB.
    """
    prefix_b = comment_prefix.encode("ascii") if comment_prefix is not None else None
    with try_open_bin(path, "rb") as fh:
        for line in fh:
            if prefix_b is not None and line.startswith(prefix_b):
                continue
            if line.endswith(b"\n"):
                line = line[:-1]
                if line.endswith(b"\r"):
                    line = line[:-1]
            if not line:
                continue
            yield line.split(b"\t")


def iter_tab_rows_str(path, comment_prefix=None, encoding="utf-8", errors="replace"):
    """
    Yield rows as list[str] split on TAB.
    """
    for parts in iter_tab_rows_bytes(path, comment_prefix=comment_prefix):
        yield [p.decode(encoding, errors) for p in parts]


def _b2s(x, encoding="utf-8", errors="replace"):
    return x.decode(encoding, errors) if isinstance(x, (bytes, bytearray)) else x


def describe(text, width=80, margin=2):
    margin = " " * margin
    # remove flanking white space
    text = text.strip()
    # insert script name 
    text = text.format(SCRIPT=os.path.split(sys.argv[0])[1])
    lines = text.split("\n")
    # title
    rule = "=" * width
    newlines = [rule, margin + lines[0], rule, "\n"]
    newline = margin
    # description
    for line in lines[2:]:
        line = line.strip()
        if line == "":
            newlines.append(newline)
            newlines.append("\n")
            newline = margin
            continue
        words = line.split()
        for word in words:
            if len(word) > width:
                newlines.append(newline)
                newlines.append(word)
                newline = margin
            elif len(newline + " " + word) > width:
                newlines.append(newline)
                newline = margin + word
            else:
                newline += (" " if newline != margin else "") + word
    if len(newline) > 0:
        newlines.append(newline)
    newlines += ["\n", rule]
    return "\n".join([k for k in newlines])


def read_contig_lengths(fasta):
    data = OrderedDict()
    header = None
    with try_open_bin(fasta, "rb") as fh:
        for line in fh:
            if not line:
                continue
            line = line.strip()
            if not line:
                continue
            if line.startswith(b">"):
                header = line[1:].split(None, 1)[0].decode("utf-8", "replace")
                data[header] = 0
            else:
                if header is not None:
                    data[header] += len(line)
    return data


def write_rowdict(rowdict=None, format=None, file=None,
                  delim="\t", precision=4, empty_field="--"):
    """ write rowdict line conditioned on format """
    if rowdict is None:
        print(delim.join([k.upper() for k in format]), file=file)
    elif set(rowdict) != set(format):
        for f in format:
            print(f, rowdict.get(f, None))
        die("Format mismatch.")
    else:
        items = []
        for f in format:
            if type(rowdict[f]) in [float, np.float32, np.float64]:
                rowdict[f] = "{A:.{B}f}".format(A=rowdict[f], B=precision)
            items.append(str(rowdict[f]) if rowdict[f] != "" else empty_field)
        try:
            print(delim.join(items), file=file)
        except Exception as e:
            die("Writing row failed:", e)


class Frame(object):
    """ manipulate generic tabular data with headers """

    def __init__(self, path, encoding="utf-8", errors="replace"):
        self._path = path
        self._encoding = encoding
        self._errors = errors
        self._fh = try_open_bin(path, "rb")

        header = None
        for line in self._fh:
            if not line:
                continue
            line = line.rstrip(b"\r\n")
            if line:
                header = line
                break
        if header is None:
            die("Empty table:", path)

        self.headers = [h.decode(encoding, errors) for h in header.split(b"\t")]

    def iter_rowdicts(self):
        n = len(self.headers)
        enc, err = self._encoding, self._errors
        for line in self._fh:
            if not line:
                continue
            line = line.rstrip(b"\r\n")
            if not line:
                continue
            parts = line.split(b"\t")
            if len(parts) < n:
                parts += [b""] * (n - len(parts))
            elif len(parts) > n:
                parts = parts[:n]
            row = [p.decode(enc, err) for p in parts]
            yield {h: v for h, v in zip(self.headers, row)}

    def close(self):
        try:
            self._fh.close()
        except Exception:
            pass

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# WORKING WITH BLAST
# ---------------------------------------------------------------
# ---------------------------------------------------------------

c_blast_fields = [
    ["qseqid", str],
    ["sseqid", str],
    ["qlen", int],
    ["slen", int],
    ["length", int],
    ["qstart", int],
    ["qend", int],
    ["sstart", int],
    ["send", int],
    ["pident", float],
    ["positive", int],
    ["gaps", int],
    ["evalue", float],
    ["bitscore", float],
    ["sstrand", str],
]
c_blast_format_string = " ".join(["6"] + [fname for [fname, ftype] in c_blast_fields])
# blast default is 500, which is sometimes too small for long contigs
c_max_target_seqs = 10000

class Hit(object):
    """
    Processes the information from a single BLAST row.
    Optimized: accepts list[bytes] from iter_tab_rows_bytes to avoid decoding every field.

    Additional information is pulled into the hit object by parsing the subject sequence identifier.
    The expected format is:
    0 : unique gene id 
    1 : taxon name 
    2+: one or more functional annotations in "system=id" format
    """

    __slots__ = (
        "qseqid","sseqid","qlen","slen","length","qstart","qend","sstart","send",
        "pident","positive","gaps","evalue","bitscore","sstrand",
        "scov","qcov","ltrim","rtrim","scov_modified","waafle_score",
        "geneid","taxon","annotations"
    )

    def __init__(self, blast_row):
        if len(blast_row) != len(c_blast_fields):
            die("inconsistent blast row: {}".format(str(blast_row)))

        # strings (decode only these)
        self.qseqid  = _b2s(blast_row[0])
        self.sseqid  = _b2s(blast_row[1])
        self.sstrand = _b2s(blast_row[14])
        self.sstrand = "-" if self.sstrand == "minus" else "+"

        # numerics (parse directly from bytes/str)
        self.qlen     = int(blast_row[2])
        self.slen     = int(blast_row[3])
        self.length   = int(blast_row[4])
        self.qstart   = int(blast_row[5])
        self.qend     = int(blast_row[6])
        self.sstart   = int(blast_row[7])
        self.send     = int(blast_row[8])
        self.pident   = float(blast_row[9])
        self.positive = int(blast_row[10])
        self.gaps     = int(blast_row[11])
        self.evalue   = float(blast_row[12])
        self.bitscore = float(blast_row[13])

        # derived coverage stats
        self.scov = (abs(self.send - self.sstart) + 1) / float(self.slen)
        self.qcov = (abs(self.qend - self.qstart) + 1) / float(self.qlen)

        # special scoverage that doesn't penalize hanging off contig end
        if self.sstrand == "-":
            sstart = self.slen - self.sstart + 1
            send   = self.slen - self.send + 1
        else:
            sstart = self.sstart
            send   = self.send

        self.ltrim = max(0, sstart - self.qstart)
        self.rtrim = max(0, self.slen - sstart - self.qlen + self.qstart)

        denom = float(self.slen - self.ltrim - self.rtrim)
        #if denom <= 0:
        #    denom = 1.0
        self.scov_modified = (send - sstart + 1) / denom
        # waafle score
        self.waafle_score = self.scov_modified * self.pident / 100.0

        # values extracted from subject header
        items = self.sseqid.split("|")
        if len(items) < 2:
            die("bad subject id header:", self.sseqid)
        self.geneid = items[0]
        self.taxon  = items[1]

        # process possible annotations
        self.annotations = {}
        if len(items) > 2:
            for k in items[2:]:
                system, name = k.split("=")
                self.annotations[system] = name


def iter_hits(blastoutfile):
    """
    Iterate through the hits in a blast file
    """
    for row in iter_tab_rows_bytes(blastoutfile, comment_prefix=None):
        yield Hit(row)


def iter_contig_hits(blastoutfile):
    """
    Iterate through hits by contig (assumes file is sorted by query)
    """
    contig, hits = None, []
    for row in iter_tab_rows_bytes(blastoutfile, comment_prefix=None):
        hit = Hit(row)
        if contig is not None and hit.qseqid != contig:
            yield contig, hits
            # reset
            hits = []
        contig = hit.qseqid
        hits.append(hit)
    # last case clean up
    if contig is not None:
        yield contig, hits


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# WORKING WITH GFFS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_gff_fields = [
    ["seqname", str],
    ["source", str],
    ["feature", str],
    ["start", int],
    ["end", int],
    ["score", float],
    ["strand", str],
    ["frame", str],
    ["attribute", str],
]

# ---------------------------------------------------------------
# gff line object (locus)
# ---------------------------------------------------------------

class Locus(object):

    def __init__(self, gff_row, attach_annotations=True):
        # gff fields
        if len(gff_row) != len(c_gff_fields):
            die("Bad GFF row:", gff_row)
        for (fname, ftype), value in zip(c_gff_fields, gff_row):
            setattr(self, fname, ftype(value) if value != "." else value)
        # annotations
        self.annotations = {}
        self.annotation_scores = {}
        for item in self.attribute.split("; "):
            match = re.search(r"^(.*?) \"(.*)\"$", item)
            if attach_annotations and match:
                system, value = match.groups()
                self.annotations[system] = value
                self.annotation_scores[system] = None
        # no name by default
        self.name = None
        self.code = ":".join([str(self.start), str(self.end), self.strand])
        #ignore status
        self.ignore = False

    def __len__(self):
        return abs(self.end - self.start) + 1

# ---------------------------------------------------------------
# gff utils
# ---------------------------------------------------------------

"""
Note: Currently we'll ignore any annotations in a user-supplied
GFF and focus on those assigned by WAAFLE.
Could update this later.
"""

def iter_loci(gff_file, attach_annotations=True):
    for row in iter_tab_rows_str(gff_file, comment_prefix="#"):
        yield Locus(row, attach_annotations=attach_annotations)


def iter_contig_loci(gff_file, attach_annotations=True):
    contig, loci = None, []
    for row in iter_tab_rows_str(gff_file, comment_prefix="#"):
        l = Locus(row, attach_annotations=attach_annotations)
        if contig is not None and l.seqname != contig:
            yield contig, loci
            loci = []
        contig = l.seqname
        loci.append(l)
    if contig is not None:
        yield contig, loci


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# WORKING WITH TAXONOMY
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_unknown = "Unknown"
c_root    = "r__Root"

# ---------------------------------------------------------------
# taxonomy object
# ---------------------------------------------------------------

class Taxonomy(object):

    def __init__(self, path):
        self.parents = {}
        self.children = {}
        for row in iter_tab_rows_str(path, comment_prefix=None):
            if len(row) < 2:
                continue
            clade, parent = row[0], row[1]
            self.parents[clade] = parent
            self.children.setdefault(parent, set()).add(clade)
        self.known_leaf_count = {}

    def get_parent(self, clade):
        return self.parents.get(clade, c_root)

    def get_children(self, clade):
        return self.children.get(clade, set())

    def get_lineage(self, clade):
        l = [clade]
        while l[-1] != c_root:
            parent = self.get_parent(clade)
            l.append(parent)
            clade = parent
        l.reverse()
        return l

    def get_lca(self, *clades):
        lca = c_root
        lineages = [self.get_lineage(c) for c in clades]
        min_depth = min([len(l) for l in lineages]) if lineages else 0
        for i in range(0, min_depth):
            level = {l[i] for l in lineages}
            if len(level) == 1:
                lca = list(level)[0]
            else:
                break
        return lca

    def get_tails(self, clades, lca):
        tails = []
        for c in clades:
            l = self.get_lineage(c)
            l.reverse()
            t = []
            for c2 in l:
                if c2 == lca:
                    break
                else:
                    t.append(c2)
            t.reverse()
            tails.append(t)
        return tails

    def get_sisters(self, clade):
        ret = set()
        parent = self.get_parent(clade)
        for child in self.get_children(parent):
            if child != clade:
                ret.add(child)
        return ret

    def get_leaf_count(self, clade):
        if clade in self.known_leaf_count:
            return self.known_leaf_count[clade]
        if clade not in self.children:
            ret = 1
        else:
            ret = 0
            for c in self.children[clade]:
                ret += self.get_leaf_count(c)
        self.known_leaf_count[clade] = ret
        return ret


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# WORKING WITH INTERVALS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

class INode(object):
    """interval node: represents an interval + some network properties"""

    __slots__ = ("start", "stop", "strand", "neighbors", "visited")

    def __init__(self, start, stop, strand="+"):
        self.start, self.stop = sorted([start, stop])
        self.strand = strand
        self.neighbors = set()
        self.visited = False

    def __len__(self):
        return self.stop - self.start + 1

    def attach(self, node):
        self.neighbors.add(node)

    def get_connected_component(self):
        # modified to use breadth-first search
        cc, front = {self}, {self}
        while any([not inode.visited for inode in front]):
            new_front = set()
            for inode in front:
                if not inode.visited:
                    inode.visited = True
                    new_front.update(inode.neighbors)
            cc.update(new_front)
            front = new_front
        return list(cc)

    def to_list(self):
        return [self.start, self.stop, self.strand]


def calc_overlap(a1, a2, b1, b2, normalize=True):
    """ compute overlap between two intervals """
    overlap = None
    a1, a2 = sorted([a1, a2])
    b1, b2 = sorted([b1, b2])
    if b1 > a2 or a1 > b2:
        overlap = 0
    else:
        outleft, inleft, inright, outright = sorted([a1, a2, b1, b2])
        overlap = inright - inleft + 1
        if normalize:
            denom = min((a2 - a1 + 1), (b2 - b1 + 1))
            overlap /= float(denom)
    return overlap


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# WORKING WITH SAM OUTPUT
# ---------------------------------------------------------------
# ---------------------------------------------------------------

"""
@HD  VN:1.0                SO:unsorted
@SQ  SN:SRS011061_k119_3   LN:961
@SQ  SN:SRS011061_k119_5   LN:837
@SQ  SN:SRS011061_k119_11  LN:502
...
61NLYAAXX100508:5:100:10002:9010   83   SRS011061_k119_37319  25069  42  101M      =  25052  -118
61NLYAAXX100508:5:100:10002:9010   163  SRS011061_k119_37319  25052  42  46M1I43M  =  25069  118
61NLYAAXX100508:5:100:10002:9075   83   SRS011061_k119_14610  17113  42  101M      =  16942  -272
61NLYAAXX100508:5:100:10002:9075   163  SRS011061_k119_14610  16942  42  46M1I43M  =  17113  272
61NLYAAXX100508:5:100:10003:17250  77   *                     0      0   *         *  0      0
61NLYAAXX100508:5:100:10003:17250  141  *                     0      0   *         *  0      0
61NLYAAXX100508:5:100:10003:3146   99   SRS011061_k119_83764  304    42  101M      =  366    152
61NLYAAXX100508:5:100:10003:3146   147  SRS011061_k119_83764  366    42  90M       =  304    -152
"""

class SAMHit(object):
    """ Some data about an aligned read in a SAM file """

    __slots__ = ("qseqid", "sseqid", "sstart", "send")

    def __init__(self, qseqid, sseqid, sstart, send):
        self.qseqid = qseqid
        self.sseqid = sseqid
        self.sstart = sstart
        self.send   = send


def cigar_length(cigar):
    counts = [int(c) for c in re.split("[A-Z]+", cigar) if c != ""]
    sigils = [s for s in re.split("[0-9]+", cigar) if s != ""]
    # ignore read-only bands
    return sum([c for c, s in zip(counts, sigils) if s in "DHMNSX="])


def iter_sam_hits(sam_file, encoding="utf-8", errors="replace"):
    """
    Very fast SAM iterator for .sam or .sam.gz:
    - reads bytes
    - skips header lines '@'
    - splits only required fields (0,2,3,5)
    """
    with try_open_bin(sam_file, "rb") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(b"@"):
                continue
            line = line.rstrip(b"\r\n")
            if not line:
                continue

            # qname, flag, rname, pos, mapq, cigar, ...
            parts = line.split(b"\t", 6)
            if len(parts) < 6:
                continue

            rname = parts[2]
            if rname == b"*":
                continue

            try:
                qname_s = parts[0].decode(encoding, errors)
                rname_s = rname.decode(encoding, errors)
                sstart  = int(parts[3])
                cigar   = parts[5].decode("ascii", "replace")
                send    = sstart + cigar_length(cigar) - 1
            except Exception:
                continue

            yield SAMHit(qname_s, rname_s, sstart, send)


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# TESTING
# ---------------------------------------------------------------
# ---------------------------------------------------------------

if __name__ == "__main__":
    pass
