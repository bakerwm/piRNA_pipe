#!/usr/bin/env python

"""
Parse the i7 and inline barcode from sequencing reads

support: TruSeq library (64nt)
# p7a = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' # 34nt 
# p7b = 'ATCTCGTATGCCGTCTTCTGCTTG' # 24nt
# p7 # 6nt



output format:
1. P7 structure: P7a-i7-P7b
"group", "total", "count", "pct", "sample"
#P7a_i7 17107874        16313570        95.4    RNAseq_shWhite_0-2h_rep1_1

2. i7 and bc pct: 
seq  i7_rc   i7_name   demx  demx_pct group total i7 i7_pct sample
>CAACTA CAACTA  TruSeq_Index29  16210525         99.4   i7      17107874        16313570        95.4    RNAseq_shWhite_0-2h_rep1_1

"""


import os
import sys
import re
import pyfastx
from hiseq.demx.sample_sheet import HiSeqIndex # check index name/seq
from collections import Counter

def top_seq(x, n=3, revcmp=False):
    t = len(x)
    a = Counter(x).most_common(n)
    out = []
    for i in a:
        k, v = i
        # format:
        rp = rev_comp(k) if revcmp else k
        rn = HiSeqIndex(rp).name # NULL if not found
        if rn is None or rn == 'NULL':
            rn = rp
        # seq, seq-rev, seq-name, pct, count
        out.append('{}\t{}\t{}\t{}\t{:5.1f}'.format(k, rp, rn, v, v/t*100))
    return out


def rev_comp(x):
    """Reverse complement DNAseq"""
    t = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A',
        'N': 'N',
    }
    s = x[::-1] # reverse
    return ''.join([t.get(i, i) for i in s]) # complement


## patterns
def parse_i7(fx, outdir, save_seq=False):
    fname = os.path.basename(fx)
    fname = fname.replace('.fq.gz', '')
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    ## output
    out_bc = os.path.join(outdir, fname + '.barcode.txt')
    out_msg = os.path.join(outdir, fname + '.stat.log')
    ## adapters
    # p1 = re.compile('AGATCG[ACGTN]{22}AGTCAC')  # p7a
    p1 = re.compile('([ACGTN]{6})AGATCG[ACGTN]{22}AGTCAC([ACGTN]{6})(.*)') # p7a+p7 
    p2 = re.compile('AGATCG[ACGTN]{22}AGTCAC([ACGTN]{6})ATCTCG[ACGT]{12}TGCTTG') # p7a+p7+p7b
    p3 = re.compile('([ACGTN]{6})ATCTCG[ACGT]{12}TGCTTG') # p7+p7b
    ## count reads
    n0 = 0
    n1 = 0
    n2 = 0
    n3 = 0
    bc = []
    i7 = []
    with open(out_bc, 'wt') as w:
        for _,s,_,_ in pyfastx.Fastx(fx):
            t = []
            n0 += 1
            if n0%1e6 == 0:
                print('{} processed'.format(n0))            
            s1 = p1.search(s)
            s3 = p3.search(s)
            if s1:
                bc.append(s1.group(1))
                i7.append(s1.group(2))
                t.extend([s1.group(1), s1.group(2), s1.group(3)])
                n1 += 1
                s2 = p2.search(s)        
                if s2:
                    n2 += 1
                    n3 += 1
            elif s3:
                t.append(s3.group(1))
                n3 += 1
            if len(t) > 0 and save_seq:
                w.write('\t'.join(t)+'\n')
    ## top i7
    kbc = '\n'.join(['>{}\tbc\t{}\t{}\t{:.1f}\t{}'.format(i, n0, n1, n1/n0*100, fname) for i in top_seq(bc, 3, revcmp=True)])
    ki7 = '\n'.join(['>{}\ti7\t{}\t{}\t{:.1f}\t{}'.format(i, n0, n1, n1/n0*100, fname) for i in top_seq(i7, 3)])
    ## message
    msg = '\n'.join([
        '-'*80,
        '{:>11}: {:>12}'.format('fastq', fx),
        '{:>11}: {:12,}'.format('total', n0),
        '\t'.join(["group", "total", "count", "pct", "sample"]),
        '#{}\t{}\t{}\t{:.1f}\t{}'.format('P7a_i7', n0, n1, n1/n0*100, fname),
        '#{}\t{}\t{}\t{:.1f}\t{}'.format('P7a_7b', n0, n2, n2/n0*100, fname),
        '#{}\t{}\t{}\t{:.1f}\t{}'.format('i7_P7b', n0, n3, n3/n0*100, fname),
        '-'*80,
        '\t'.join(['demx', 'demx_rc', 'demx_name', 'demx_n', 'demx_pct', 'group', 'total', 'i7', 'i7_pct', 'sample']),
        ki7,
        kbc,
#         '[{:<}]: \n{}'.format('i7', ki7),
#         '[{:<}]: \n{}'.format('barcode', kbc),
        '-'*80,
    ])
    print(msg)
    with open(out_msg, 'wt') as w:
        w.write(msg+'\n')


def main():
    if len(sys.argv) < 3:
        print('Usage: python parse_i7.py <outdir*> <in.fq*> <write_seq: 1|0>')
        sys.exit(1)
    outdir, fx = sys.argv[1:3]
    save_seq = len(sys.argv) == 4
    parse_i7(fx, outdir, save_seq)
    
    
if __name__ == '__main__':
    main()
        
# 