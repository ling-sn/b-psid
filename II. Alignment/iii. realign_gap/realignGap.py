#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Created: 2021-11-03 16:29
"""
Python script that realigns .bam files.
"""
from pathlib import Path
import traceback
import argparse
import re
import parasail
import pysam
import subprocess

## Parasail and cigar functions
re_split_cigar = re.compile(r"(?P<len>\d+)(?P<op>\D+)")

def sam2align(record):
    query_string = ""
    match_string = ""
    ref_string = ""
    for query_pos, ref_pos, ref_base in record.get_aligned_pairs(with_seq=True):
        if ref_base is None and query_pos is None:
            query_string += "."
            match_string += " "
            ref_string += "N"
        elif query_pos is None and ref_pos is not None:
            query_string += "-"
            match_string += " "
            ref_string += ref_base
        elif ref_pos is None and query_pos is not None:
            query_base = record.query_sequence[query_pos]
            query_string += query_base
            match_string += " "
            ref_string += "-"
        elif ref_pos is not None and query_pos is not None:
            query_base = record.query_sequence[query_pos]
            query_string += query_base
            if query_base == ref_base:
                match_string += "|"
            else:
                match_string += " "
            ref_string += ref_base
        else:
            pass
    return "\n".join([query_string, match_string, ref_string])

def cigar_ops_from_start(cigar):
    """
    Yield cigar string operations from start of cigar (in order).

    :param cigar: cigar string.
    :yields: str, op. length, str op. type
    """
    for m in re.finditer(re_split_cigar, cigar):
        yield m.group("len"), m.group("op")

def parasail_to_sam(result, seq):
    """
    Extract reference start and sam compatible cigar string.

    :param result: parasail alignment result.
    :param seq: query sequence.
    :returns: reference start coordinate, cigar string.
    """
    cigstr = result.cigar.decode.decode()

    first = next(cigar_ops_from_start(cigstr))
    prefix = "".join(first)
    rstart = result.cigar.beg_ref
    cliplen = result.cigar.beg_query
    clip = "" if cliplen == 0 else "{}S".format(cliplen)
    if first[1] == "I":
        pre = "{}S".format(int(first[0]) + cliplen)
    elif first[1] == "D":
        pre = clip
        rstart = int(first[0])
    else:
        pre = "{}{}".format(clip, prefix)

    mid = cigstr[len(prefix) :]
    end_clip = len(seq) - result.end_query - 1
    suf = "{}S".format(end_clip) if end_clip > 0 else ""
    new_cigstr = "".join((pre, mid, suf))
    return rstart, new_cigstr

def parasail_alignment(query, ref):
    """
    Run a Smith-Waterman alignment between two sequences.

    :param query: the query sequence.
    :param ref: the reference sequence.
    :returns: reference start co-ordinate, cigar string
    """
    result = parasail.sw_trace_striped_32(query, ref, 3, 2, parasail.dnafull)
    rstart, cigar = parasail_to_sam(result, query)
    return rstart, cigar

def get_splicing_site_from_cigar(cigar):
    """
    Get splicing site from cigar string.

    :param cigar: cigar string.
    :returns: splicing site.
    """
    splicing_site = []
    n = 0
    for x, y in cigar_ops_from_start(cigar):
        x = int(x)
        if y == "N":
            # if n > 0:
            splicing_site.append((n, x))
            # n = 0
        elif y in "X=MDS":
            n += x

    return splicing_site

def degenerate_cigar(cigar):
    new_cigar = ""
    n = 0
    for x, y in cigar_ops_from_start(cigar):
        x = int(x)
        if y not in "X=":
            if n > 0:
                new_cigar += f"{n}M"
            new_cigar += f"{x}{y}"
            n = 0
        else:
            n += x
    if n > 0:
        new_cigar += f"{n}M"
    return new_cigar

def split_cigar(cigar, split_points):
    """
    Input the cigar string before degenration.

    52= for example.
    """
    new_cigar = ""
    split_points.append((float("inf"), 0))
    w = 0
    p = 0
    s, l = split_points.pop(0)
    s2, _ = split_points[0]

    for x, y in cigar_ops_from_start(cigar):
        x = int(x)
        if y in "X=MS":
            p += x
        if p <= s:
            new_cigar += f"{x}{y}"
            w += x
        else:
            while x > 0:
                s2, _ = split_points[0]
                if p - s <= s2 - s:
                    new_cigar += f"{x - p + s}{y}{l}N{p - s}{y}"
                    x -= s2
                else:
                    x -= s - w
                    new_cigar += f"{x}{y}{l}N"
                    w = 0
                s, l = split_points.pop(0)
    return new_cigar

def run_realign(input_bam_name, output_bam_name, fasta_dir, discard):
    fafile = pysam.FastaFile(fasta_dir) ## specify input FASTA file
    bamfile = pysam.AlignmentFile(input_bam_name, "rb")
    outfile = pysam.AlignmentFile(output_bam_name, "wb", template=bamfile)
    
    try:
        for read in bamfile.fetch(until_eof=True):
            if (
                "D" in read.cigarstring
                or "S" in read.cigarstring
                or not read.get_tag("MD").isnumeric()
            ) and read.reference_name in fafile.references:
                if "N" not in read.cigarstring:
                    align_start = max(read.reference_start - 20, 0)
                    align_end = read.reference_end + 20
                    ref = fafile.fetch(read.reference_name, align_start, align_end)
                    s, c = parasail_alignment(read.query_sequence, ref)
                    # update align position
                    read.reference_start = align_start + s
                    # use a old cigar spec...
                    c2 = degenerate_cigar(c)
                    if read.cigarstring != c2:
                        # write original cigar
                        origin_cigar = read.cigarstring
                        read.set_tag("OC", origin_cigar)
                        # update cigar string
                        read.cigarstring = c2
                # realign splicing reads
                else:
                    ref = ""
                    align_start = max(read.reference_start - 20, 0)
                    align_end = read.reference_end + 20
                    ref += fafile.fetch(
                        read.reference_name, align_start, read.reference_start
                    )
                    ref += read.get_reference_sequence()
                    ref += fafile.fetch(
                        read.reference_name, read.reference_end, align_end
                    )
                    s, c = parasail_alignment(read.query_sequence, ref)
                    # update align position
                    read.reference_start = align_start + s
                    split_points = get_splicing_site_from_cigar(read.cigarstring)
                    c2 = degenerate_cigar(split_cigar(c, split_points))
                    if read.cigarstring != c2:
                        # write original cigar
                        origin_cigar = read.cigarstring
                        read.set_tag("OC", origin_cigar)
                        # update cigar string
                        read.cigarstring = c2
            if read.infer_query_length() != read.query_length:
                if not discard:
                    # revert origin
                    read.reference_start = align_start
                    read.cigarstring = origin_cigar
                    read.set_tag("OC", None)
                    outfile.write(read)
            else:
                outfile.write(read)

    except Exception as e:
        print(f"Failed to realign {input_bam_name.name}: {e}")
        traceback.print_exc()
        raise

def bam_index(output_bam_name):
    sorted_bam = output_bam_name.with_suffix(".sorted.bam")
    try:
        subprocess.run(["samtools", "sort", 
                        "-o", str(sorted_bam),
                        str(output_bam_name)],
                        check = True,
                        capture_output = True,
                        text = True)
        subprocess.run(["samtools", "index", str(sorted_bam)], ## create .bai from .bam
                        check = True,
                        capture_output = True,
                        text = True)
        subprocess.run(["rm", str(output_bam_name)],
                        check = True,
                        capture_output = True,
                        text = True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to convert {output_bam_name.name} to .bai: {e}")
        print("STDERR:", e.stderr)
        print("STDOUT:", e.stdout)
        traceback.print_exc()
        raise

## Main function
def main(star_folder, subf, fasta_dir, discard):
    current_path = Path.cwd()
    subfolder = current_path/star_folder/subf
    
    try: 
        if subfolder.is_dir():
            processed_folder = current_path/"realignments"/f"{subfolder.stem}"
            processed_folder.mkdir(exist_ok = True, parents = True)
            
            ## Obtain dedup .bam
            input_bam_name = next(bam for bam in subfolder.glob("*_dedup.bam"))
            output_bam_name = processed_folder/f"{input_bam_name.stem}_realigned.bam"

            ## Sort .bam if not indexed already
            bai_exists = list(subfolder.glob("*.bai"))
            if not bai_exists:
               subprocess.run(["samtools", "index", str(input_bam_name)], ## create .bai from .bam
                               check = True,
                               capture_output = True,
                               text = True)

            ## Run realignment
            run_realign(input_bam_name, output_bam_name, fasta_dir, discard)
            bam_index(output_bam_name)

    except Exception as e:
        print(f"Failed to obtain BAM files for realignment: {e}")
        traceback.print_exc()
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Realigns BAM files in BID-Seq pipeline.")
    parser.add_argument("--star_folder", help = "Name of overall STAR aligned folder containing subfolders",
                        required = True)
    parser.add_argument("--subf", help = "Name of STAR aligned subfolder that you wish to realign", 
                        required = True)
    parser.add_argument("--fasta_dir", help = "Path of reference fasta file", required = True)
    parser.add_argument("--discard", action = argparse.BooleanOptionalAction, default = False, 
                        help = "Write discarded reads into file for debugging. Use '--no-discard' to disable and '--discard' to enable.")
    args = parser.parse_args()

    print("Starting realignment...")
    main(args.star_folder, args.subf, args.fasta_dir, args.discard)
    print("Realignment finished.")