#!/usr/bin/env python

"""
Convert bam to bed

save query-sequence in name:
name:seq
"""

import os
import sys
import pysam 
import fileinput


def bam2bed(fh):
    for read in fh:
        if read.is_unmapped or \
            read.is_supplementary or \
            read.is_secondary:
            next
        # output
        bed = [
            read.reference_name,
            read.reference_start,
            read.reference_end,
            '{}:{}'.format(read.query_name, read.query_sequence),
            '255',
            '-' if read.is_reverse else '+'
        ]
        bed = list(map(str, bed)) # to str
        yield bed
        

def main():
    if len(sys.argv) < 2:
        print('Usage: python bam2bed.py <bam>')
        sys.exit(1)
    s = sys.argv[1]
    # check bam and bai
    if not s.endswith('.bam'):
        print('not a bam file: {}'.format(s))
        sys.exit(1)
    if not os.path.exists(s+'.bai'):
        pysam.index(s)
#     # read data
#     s = sys.argv[1] if len(sys.argv) > 1 else '-' # stdin
    sam = pysam.AlignmentFile(s, "r")
    # output
    for b in bam2bed(sam):
        print('\t'.join(b))
    
    
if __name__ == '__main__':
    main()

#