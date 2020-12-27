#!/usr/bin/env python3

"""
Extract peaks from TE

Check two criteria: peaks_length, peak_reads

input: bam (convert to bedGraph)
output: peak, stat

"""

import os
import sys
import pandas as pd
import numpy as np
import pybedtools


def bed2fa(bed, ref=None, output=None):
    """
    Extract fasta sequence from reference, using BED
    """
    if isinstance(ref, str) and isinstance(output, str):
        pybedtools.BedTool(bed).sequence(fi=ref, fo=output, name=True, tab=True)
    else:
        print('ref=, output=, required')


def bam2bdg(bam, bdg, strand='+'):
    """
    Convert bam to bdg
    strand = +, -

    a = pybedtools.BedTool('a.bam')
    b = a.genome_coverage(bg=True, strand='+')
    """
    if os.path.exists(bdg):
        print('file exists, bam2bdg skipped: {}'.format(bdg))
    else:
        a = pybedtools.BedTool(bam)
        b = a.genome_coverage(bg=True, strand=strand).saveas(bdg)


def bdg2peak(bdg, peak, cutoff=5, strand='.'):
    """
    Use macs2 bdgpeakcall to extrack peaks
    """
    # save peaks (unstranded)
    peak_tmp = os.path.splitext(peak)[0] + '.unstranded' + '.bed'
    cmd = 'macs2 bdgpeakcall -i {} -c {} -l 23 -g 10 -o {}'.format(bdg, cutoff, peak_tmp)

    if os.path.exists(peak_tmp):
        print('file exists, bdg2peak skipped, {}'.format(peak_tmp))
    else:
        os.system(cmd)

    # Add strand to peak
    if os.path.exists(peak):
        print('file exists, bdg2peak fix, skipped, {}'.format(peak))
    else:
        p1 = pybedtools.BedTool(peak_tmp).to_dataframe().drop(0) # remove 1st row
        p1.strand = strand # fix
        # fix peak name
        p1.name = [i.replace(peak_tmp, '') for i in p1.name.tolist()]
        p1.name = p1.chrom + p1.name
        # fix type
        p1 = p1.astype({'start': 'int', 'end': 'int', 'thickStart': 'int', 'thickEnd': 'int'})
        p2 = pybedtools.BedTool.from_dataframe(p1).saveas(peak, trackline=None)


def peak2count(peak, bam, output):
    """
    Count reads for each peak
    intersect

    a = pybedtools.BedTool('a.peak')
    b = pybedtools.BedTool('a.bam')
    a_and_b = a.intersect(b, c=True) # last column
    """
    if os.path.exists(output):
        print('file exists, peak2count skipped: {}'.format(output))
    else:
        a = pybedtools.BedTool(peak) # bed
        b = pybedtools.BedTool(bam)
        ab = a.intersect(b, c=True, s=True).saveas(output)


def peak2stat(peak, by_te=False):
    """
    Count peak reads, peaks, length
    return: reads, peak_length
    """
    if isinstance(peak, pd.DataFrame):
        df = peak
    elif isinstance(peak, str):
        df = pybedtools.BedTool(peak).to_dataframe() 
    else:
        df = pd.DataFrame(columns = ['chrom', 'start', 'end', 'blockSizes'])

    # summarize
    df['length'] = df.end - df.start #

    if by_te:
        df = df.groupby('chrom')

    npeak = df.count().start # number of peaks
    nread = df.sum().blockSizes # number of reads
    nlen = df.sum().length

    return [npeak, nread, nlen]


def get_cutoff(peak, num=10):
    """
    auto get the cutoff based on the bdg file
    from = quantile(0.5)
    to = quantile(0.9)
    num = 4
    """
    a = pybedtools.BedTool(peak).to_dataframe() # chrom, start, end, name(count)
    start = a.name.quantile(0.5)
    end = a.name.quantile(0.99)
    return np.linspace(start, end, num).astype(int).tolist()


def scan_peak(bam, outdir, strand='+', num=10, ref=None):
    """
    scan peaks on TEs
    return the peak files, stat file
    """
    # filename
    fname = os.path.splitext(os.path.basename(bam))[0] # name
    suffix = 'rev' if strand == '-' else 'fwd'

    # output
    outdir = os.path.join(outdir, suffix) # subdir_1
    if not os.path.exists(outdir):
            os.makedirs(outdir) # make dir

    # convert bam to bedGraph
    bdg = os.path.join(outdir, fname + '.bdg')
    bam2bdg(bam, bdg, strand=strand)

    # scan cutoff
    dn = []
    dn_te = []
    cutoff_list = get_cutoff(bdg, num)
    for cutoff in cutoff_list:
        outdir_sub = os.path.join(outdir, 'cutoff_{}'.format(cutoff))
        if not os.path.exists(outdir_sub):
            os.makedirs(outdir_sub) # make dir
        peak = os.path.join(outdir_sub, fname + '.peak.bed')
        peak2 = os.path.join(outdir_sub, fname + '.peak.count.bed')
        peak_fa = os.path.join(outdir_sub, fname + '.peak.fa')
        bdg2peak(bdg, peak, cutoff, strand) # strand fixed
        bed2fa(peak, ref, peak_fa) # convert to fasta
        peak2count(peak, bam, peak2)
        # stat for all TEs
        npeak, nread, nlen = peak2stat(peak2)
        dn.append([npeak, nread, nlen, cutoff])
        # stat for each TE
        npeak2, nread2, nlen2 = peak2stat(peak2, True)
        df_n = pd.concat([npeak2, nread2, nlen2], axis = 1).reset_index()
        df_n.columns = ['chrom', 'npeak', 'nread', 'nlen']
        df_n['cutoff'] = cutoff
        dn_te.append(df_n)

    # convert to dataframe, and save to file
    stat = os.path.join(outdir, 'cutoff_scan.stat')
    df = pd.DataFrame(dn, columns = ['npeak', 'nread', 'nlen', 'cutoff']).to_csv(stat, index=False)

    # convert to dataframe, 
    stat_te = os.path.join(outdir, 'cutoff_scan.te.stat')
    df = pd.concat(dn_te, axis = 0).to_csv(stat_te, index=False)


def main():
    if len(sys.argv) < 4:
        sys.exit('get_peaks.py in.bam ref.fa outdir')

    bam = sys.argv[1]
    ref = sys.argv[2]
    outdir = sys.argv[3]

    scan_peak(bam, outdir, strand='+', ref=ref)
    scan_peak(bam, outdir, strand='-', ref=ref)


if __name__ == '__main__':
    main()


# EOF