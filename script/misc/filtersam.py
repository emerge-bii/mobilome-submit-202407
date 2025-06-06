#!/usr/bin/env python3
# coding: utf-8

"""
Tools to filter SAM/BAM files by percent identity and percent of matched sequence
"""

import sys
import os
import argparse
import re
import subprocess
from pathlib import Path

import pysam
from parallelbam.parallelbam import parallelizeBAMoperation

def sumMatchesAndMismatches(segment):
    """
    Get total matches/mismatches from CIGAR string (M field)
    Code dictionary:
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9
    """
    return sum(
        [value for (code, value) in segment.cigartuples if code in {0, 7, 8}]
    )

def getNumberOfMatches(segment):
    """
    Get numnber of matches from alignment
    Do not consider insertion/deletion as mismatches
    """
    parsed_MD = segment.get_aligned_pairs(with_seq=True)
    return len([
        base for (read_pos, ref_pos, base) in parsed_MD 
        if ((base is not None and base.isupper()) and read_pos is not None)
    ])

def getQueryLength(segment):
    """
    Compute query length from CIGAR field corresponding to
    query sequence. 
    Following: https://stackoverflow.com/questions/39710796/infer-the-length-of-a-sequence-using-the-cigar
    
    NOTE2SELF: check difference w/ 'pysam.AlignedSegment.infer_query_length'
    
    Cigar fields which 'consume sequence': M, I, S, =, X
    
    Code dictionary:
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9
    """
    codes = [0, 1, 4, 7, 8]
    return sum(
        [value for (code, value) in segment.cigartuples if code in codes]
    )

def percent_matched(segment):
    """
    Compute percentage of sequence that has been matched to reference
    """
    seq_length = getQueryLength(segment)
    n_matches = getNumberOfMatches(segment)
    return 100 * (n_matches / seq_length)
    
def percent_identity(segment):
    """
    Compute percent identity from MD tag of aligned segment.
    segment: pysam AlignedSegment object.
    """
    try:
        identity = 100 * (getNumberOfMatches(segment) / sumMatchesAndMismatches(segment)) 
    except ZeroDivisionError as e:
        identity = 0
    return identity

def percent_aligned(segment):
    """
    Compute percentage of sequence that has been aligned to reference
    """
    seq_length = getQueryLength(segment)

    aligned_length = sumMatchesAndMismatches(segment)
    return 100 * (aligned_length / seq_length)

def has_MD_tag(segment):
    return 'MD' in [tag for (tag, _) in segment.get_tags()]

def filterSAMbyIdentity(input_path: Path, output_path: Path = None,
                        identity_cutoff: float = 95.0) -> None:
    """
    Filter aligned segments in BAM or SAM file with percent identity
    equal or above identity_cutoff value.
    """
    input_path = Path(input_path)
    #file_ext = re.search('.(s|b)am', input_path.as_posix()).group()
    if output_path is None:
        output_path = '-'
    output_path = Path(output_path)
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    pysam.set_verbosity(save)                                      
    for segment in samfile:
        if (has_MD_tag(segment) and percent_identity(segment) >= identity_cutoff):
            filtered_sam.write(segment)           
    filtered_sam.close()
    samfile.close()

def filterSAMbyPercentMatched(input_path: Path, output_path: Path = None,
                              matched_cutoff: float = 80.0) -> None:
    """
    Filter aligned segments in BAM or SAM file with percent of matched
    based equal or higher than matched_cutoff. 
    
    Percent of matched bases is computed as the fraction of matches in
    the total query length.
    """
    input_path = Path(input_path)
    #file_ext = re.search('.(s|b)am', input_path.as_posix()).group()
    if output_path is None:
        output_path = '-'
    output_path = Path(output_path)
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    pysam.set_verbosity(save)
    for segment in samfile:
        if (has_MD_tag(segment) and percent_matched(segment) >= matched_cutoff):
            filtered_sam.write(segment)
            
    filtered_sam.close()
    samfile.close()

def filterSAMbyPercentAligned(input_path: Path, output_path: Path = None,
                              aligned_cutoff: float = 85.0) -> None:
    """
    Filter aligned segments in BAM or SAM file with aligned
    percent
    
    Percent of matched bases is computed as the fraction of aligned in
    the total query length.
    """
    input_path = Path(input_path)
    #file_ext = re.search('.(s|b)am', input_path.as_posix()).group()
    if output_path is None:
        output_path = '-'
    output_path = Path(output_path)
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    pysam.set_verbosity(save)
    for segment in samfile:
        if (has_MD_tag(segment) and percent_aligned(segment) >= aligned_cutoff):
            filtered_sam.write(segment)
            
    filtered_sam.close()
    samfile.close()
    
def filterSAMbyPercentAlignedAndIdentity(input_path: Path, output_path: Path = None,
                              identity_cutoff: float = 95.0, aligned_cuoff: float = 85.0) -> None:
    """
    Filter aligned segments in BAM or SAM file with both percent of aligned
    and percent identity
    
    Percent of matched bases is computed as the fraction of matches in
    the total query length.
    """
    input_path = Path(input_path)
    #file_ext = re.search('.(s|b)am', input_path.as_posix()).group()
    if output_path is None:
        output_path = '-'
    output_path = Path(output_path)
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(input_path, 'r')
    filtered_sam = pysam.AlignmentFile(output_path, 'w', template=samfile)
    pysam.set_verbosity(save)
    for segment in samfile:
        if (has_MD_tag(segment) and percent_aligned(segment) >= aligned_cutoff and percent_identity(segment) >= identity_cutoff):
            filtered_sam.write(segment)
            
    filtered_sam.close()
    samfile.close()

def filterSAM(input_path: Path, output_path: Path = None,
              filter_by: str = 'identity', cutoff: float = 95.0,
              n_processes: int = None) -> None:
    """
    Filter aligned segments in BAM or SAM file by percent identity or percent
    of matched sequence.
    """
    if filter_by == 'identity':
        filter_method = filterSAMbyIdentity
    elif filter_by == 'matched':
        filter_method = filterSAMbyPercentMatched
    elif filter_by == 'aligned':
        filter_method = filterSAMbyPercentAligned
    else:
        raise ValueError('Invalid filter, available filters are: "identity", "matched", "aligned"')

    if not (cutoff >=0 and cutoff <= 100):
        raise ValueError('Cutoff value must be between 0 and 100.')
    
    if n_processes is None:
        filter_method(input_path, output_path, cutoff)
    #else:
    #    if '.sam' in input_path.name:
    #        sys.stderr.write('[INFO] Converting sam file to bam for processing\n')
    #        bam_path = Path(os.fspath(input_path).replace('.sam', '.bam'))
    #        sam2bam(input_path, bam_path)
    #    parallelizeBAMoperation(path_to_bam=bam_path, callback=filter_method,
    #                            callback_additional_args=[cutoff],
    #                            n_processes=n_processes, output_path=output_path)

def sam2bam(sam_file: Path, output_dir: Path = None) -> None:
    """
    Convert sam file to bam
    """
    if output_dir is None and os.fspath(sam_file) != '-':
        output_dir = Path(os.fspath(sam_file).replace('.sam', '.bam'))
    elif os.fspath(sam_file) == '-':
        output_dir = Path('tmp.bam')

    subprocess.run(f'samtools view {sam_file} -S -b -o {output_dir}', shell=True)

def main():
    parser = argparse.ArgumentParser(
        prog='filtersam',
        description='Tools to filter SAM/BAM files by percent identity, percent aligned or percent of matched sequence',
        epilog="  \n",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('bam', type=Path,
                        help = 'path to bam/sam file with MD tag')
    parser.add_argument('-m', '--filter-by', metavar='', type=str,
                        default='identity',
                        dest='filter_by',
                        choices = ['identity', 'aligned', 'matched'],
                        help='pick one of identity (matched percent of alignment), aligned (aligned percent of query), matched (matched percent of query)')
    parser.add_argument('-c', '--cutoff', metavar='', type=float,
                        dest='cutoff',
                        default=0.95,
                        help='minimal cutoff for metric in "--filter-by"')
    parser.add_argument('-p', '--processes', metavar='', type=int,
                        dest='processes',
                        help='number of processes for parallelization')
    parser.add_argument('-o', '--output', metavar='', type=Path,
                        default=None, dest='out',
                        help='path to output file')

    args = parser.parse_args()

    cutoff = args.cutoff
    if not (cutoff >=0 and cutoff <= 100):
        raise ValueError('Cutoff value must be between 0 and 100.')
    
    #bam = Path(args.bam).absolute()
    bam = Path(args.bam)

    if not bam.is_file() and os.fspath(bam) != '-':
        sys.stderr.write(f'[ERROR] Specified bam ({bam}) file does not exist\n')
        sys.exit()

    filterSAM(input_path=bam, output_path=args.out, filter_by=args.filter_by, cutoff=args.cutoff, n_processes=None)

if __name__ == '__main__':
    main()
