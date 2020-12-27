#!/usr/bin/env python3

"""
piRNA analysis (small RNAseq) 


version: 2020-07-28
update:
1. remove temp fastq files
2. gzip fastq files

version: 2020-07-25
update:
1. collapse reads 


date: 2020-07-23
in brief:

1. remove 3' adapters  
2. remove structural RNAs (tRNAs, rRNAs) (unique + multiple)     
3. filt by length, 23-29 nt   
4. collapse reads, only compare 1-23 nt (*)  
   collapse reads, (regular)  
   trim reads to 23nt from 3' end  

5. map to TE consensus (unique, multiple, unique + multiple)   
6. map to piRNA clusters (unique, multiple, unique + multiple)   
7. map to genome (not-TE, no-piRNAcluster)    
8. Overall, map to genome

"""

import os
import sys
import re
import gzip
import fnmatch
import binascii
import shutil
import logging
import argparse
import pybedtools
import pysam
import pandas as pd
import numpy as np
from xopen import xopen


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')


def get_args():
    """
    require arguments:
    -i fq
    -o outdir 
    --index-rRNA  
    --index-miRNA 
    --index-te  
    --index-piRNAcluster 
    --index-genome 
    --ref-te
    --ref-piRNAcluster
    ...
    """
    parser = argparse.ArgumentParser(description='piRNA pipe')
    parser.add_argument('-i', '--input', required=True,
        help='fastq file') 
    parser.add_argument('-o', '--output', required=True,
        help='directory to save results') 
    parser.add_argument('-c', '--collapse', action='store_true',
        help='collapse fastq reads')
    parser.add_argument('-g', '--genome', default='dm6',
        help='reference genome, default: [dm6]')
    parser.add_argument('-ir', '--index-rRNA', default='smRNA',
        help='bowtie_index for rRNA, default: ')
    parser.add_argument('-im', '--index-miRNA', default='hairpin',
        help='bowtie_index for miRNA (hairpin), default:')
    parser.add_argument('-ite', '--index-te', default='transposon',
        help='bowtie_index for te (TE consensus), default:') 
    parser.add_argument('-ip', '--index-piRNAcluster', default='piRNAcluster',
        help='bowtie_index for piRNA cluster, default:') 
    parser.add_argument('-ig', '--index-genome', default='genome',
        help='bowtie_index for genome, default:') 
    parser.add_argument('-te', '--te-fa', default='te.fa',
        help='TE sequence in fasta format') 
    parser.add_argument('-pi', '--piRNAcluster-fa', default='piRNAcluster.fa',
        help='piRNA cluster sequence in fasta format')

    return parser.parse_args()


def is_gz(x):
    if os.path.exists(x):
        with open(x, 'rb') as test_f:
            return binascii.hexlify(test_f.read(2)) == b'1f8b'
    else:
        if x.endswith('.gz'):
            return True
        else:
            return False


def gzip_cmd(src, dest=None, decompress=True, rm=True, overwrite=False):
    """
    Gzip Compress or Decompress files using gzip module in python
    rm, True/False, whether remove old file

    # check the src file by extension: .gz
    """
    if not file_exists(src):
        log.error('src file not exists: {}'.format(src))
        return None

    # output file
    if dest is None:
        if decompress:
            dest = src.rstrip('.gz')
        else:
            dest = src + '.gz'

    # run
    if os.path.exists(dest):
        log.warning('file exists, gzip_cmd() skipped: {}'.format(dest))
    else:
        if decompress:
            if is_gz(src):
                with xopen(src, 'rb') as r, xopen(dest, 'wb') as w:
                    shutil.copyfileobj(r, w)
            else:
                log.warning('src is not gzipped, gzip_cmd() skipped: {}'.format(src))
                shutil.copy(src, dest) # copy file
        else:
            if is_gz(src):
                log.warning('src is gzipped, gzip_cmd() skipped: {}'.format(src))
                shutil.copy(src, dest)
            else:
                with xopen(src, 'rb') as r, xopen(dest, 'wb', compresslevel=1) as w:
                    shutil.copyfileobj(r, w)

    # output
    if rm is True:
        os.remove(src)

    return dest


def listdir(path, full_name=True, recursive=False, include_dir=False):
    """
    List all the files within the path
    """
    out = []
    for root, dirs, files in os.walk(path):
        if full_name:
            dirs = [os.path.join(root, d) for d in dirs]
            files = [os.path.join(root, f) for f in files]
        out += files

        if include_dir:
            out += dirs

        if recursive is False:
            break

    return sorted(out)


def listfile(path='.', pattern='*', full_name=True, recursive=False):
    """
    Search files by the pattern, within directory
    fnmatch.fnmatch()

    pattern:

    *       matches everything
    ?       matches any single character
    [seq]   matches any character in seq
    [!seq]  matches any char not in seq

    An initial period in FILENAME is not special.
    Both FILENAME and PATTERN are first case-normalized
    if the operating system requires it.
    If you don't want this, use fnmatchcase(FILENAME, PATTERN).

    example:
    listfile('./', '*.fq')
    """
    fn_list = listdir(path, full_name, recursive, include_dir=False)
    fn_list = [f for f in fn_list if fnmatch.fnmatch(os.path.basename(f), pattern)]
    return sorted(fn_list)


def filename(x):
    """
    Extract the name of fastq file
    """
    if x is None:
        xname = None
    else:
        if x.endswith('.gz'):
            xname = os.path.basename(x).replace('.gz', '')
        else:
            xname = os.path.basename(x)

        xname = os.path.splitext(xname)[0]

        # remove ".unmap"
        xname = xname.replace('.unmap', '')

    return xname


def fx_type(x):
    """
    Check x is fastq or fasta
    reading the first 100 record
    fq, qual (str)
    fa, qual (None)
    """
    if not file_exists(x):
        log.info('file not found, fx_type() skipped: {}'.format(x))
        return None

    f_list = []
    i = 0
    with xopen(x, 'rt') as r:
        for n, s, q, m in readfq(r):
            i += 1
            if i > 100:
                break #
            fx = 'fa' if q is None else 'fq'
            f_list.append(fx)
    
    # check
    tag = list(set(f_list))

    if len(tag) == 1:
        return tag[0]
    else:
        log.error('unknown format: {}'.format(x))


def is_fastq(x):
    """
    check input file, in fastq format
    """
    return fx_type(x) == 'fq'


def is_fasta(x):
    """
    check input file, in fasta format
    """
    return fx_type(x) == 'fa'


def file_exists(x):
    """
    check input exists, file/path
    """
    if isinstance(x, str):
        tag = os.path.exists(x)
    elif isinstance(x, list):
        tag = all([file_exists(i) for i in x])
    else:
        tag = False # unknown

    return tag


def check_path(x):
    """
    check the dir, create dir
    """
    if os.path.isfile(x):
        pass
    elif not os.path.exists(x):
        os.makedirs(x)
    else:
        pass


def fx_collapse(x, outdir, trim_to=0, outfmt_fq=False, gzipped=True):
    """
    collapse fastx by seq
    
    outdir, directory saving the results
    trim_to, trim the fastx to specific length
    outfmt_fq, default [False], output in fasta format
    """
    fname = filename(x)
    ftype = 'fq' if outfmt_fq else 'fa'
    fx_out = os.path.join(outdir, '{}.{}'.format(fname, ftype))
    if gzipped:
        fx_out += '.gz'

    if file_exists(fx_out):
        log.info('file exists, fx_collapse() skipped: {}'.format(fname))
    else:
        # collapse
        d = {}
        with xopen(x, 'rt') as r:
            for n, s, q, m in readfq(r):
                d[s] = d.get(s, 0) + 1

        # sort by value (decreasing)
        i = 0
        with xopen(fx_out, 'wt') as w:
            for k, v in sorted(d.items(), key=lambda item: item[1], reverse=True):
                i += 1
                out = '>{}-{}\n{}\n'.format(i, v, k)
                if outfmt_fq:
                    out = '@{}-{}\n{}\n+\n{}\n'.format(i, v, k, 'I'*len(k))
                w.write(out)

    return fx_out


def fa2fq(fa, fq):
    """
    Convert fasta to fastq
    "I" for quality
    """
    fq_dir = os.path.dirname(fq)
    check_path(fq_dir)

    if file_exists(fq):
        print('file exists, fa2fq() skipped:')
    else:
        with xopen(fa, 'rt') as r, xopen(fq, 'wt') as w:
            for n, s, q, m in readfq(r):
                fq = '@{}\n{}\n+\n{}\n'.format(n, s, q)
                w.write(fq)


def fx_stat(fx, outdir=None):
    """
    seqkit stat fx
    """
    fname = filename(fx)
    ftype = fx_type(fx) # fa/fq
    fout = '{}.{}.stat'.format(fname, ftype)

    if outdir is None:
        outdir = os.path.dirname(fx)

    # default
    fname = filename(fx)
    check_path(outdir)

    # run
    # stat = os.path.join(outdir, fname + '.fx.stat')
    # stat = os.path.join(outdir, os.path.basename(fx) + '.stat')
    stat = os.path.join(outdir, fout)
    cmd = 'seqkit stat {} > {}'.format(fx, stat)

    if file_exists(stat):
        print('file exists, fx_stat() skipped')
    else:
        os.system(cmd)

    return stat


def readfq(fh): # this is a generator function
    """
    source: https://github.com/lh3/readfq/blob/master/readfq.py
    processing fastq file
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fh: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        [name, _, comment], seqs, last = last[1:].partition(" "), [], None
        for l in fh: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None, comment # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fh: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs), comment; # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None, comment # yield a fasta record instead
                break


def revComp(s):
    d = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    s = [d[c] for c in s]
    return ''.join(s[::-1])


def bam_to_bed(bam, bed=None):
    """
    Convert BAm to bed, using pybedtools
    """
    if bed is None:
        bed = os.path.splitext(bam)[0] + '.bed'

    if file_exists(bed):
        print('file exists, bam_to_bed() skipped: {}'.format(os.path.basename(bam)))
    else:
        pybedtools.BedTool(bam).bam_to_bed().saveas(bed)


def bam_to_fq(bam, fq=None):
    """
    Convert bam to fastq
    using samtools fastq 

    Note:
    multiple mapping reads

    see also: https://github.com/dpryan79/Answers/blob/master/bioinfoSE_2149/convert.py
    """
    fname = filename(bam)
    if fq is None:
        fq = os.path.splitext(bam)[0] + '.fq' # default

    # output fmt
    fx = fq
    if fx.endswith('.gz'):
        fx = fx.rstrip('.gz')
    ext = os.path.splitext(fx)[1]
    ftype = 'fa' if ext.lower() in ['.fa', '.fasta'] else 'fq'

    # read bam
    samfile = pysam.AlignmentFile(bam, 'rb')

    if file_exists(fq):
        print('file exists, bam_to_fq() skipped: {}'.format(os.path.basename(fq)))
    else:
        of = xopen(fq, 'wt')
        for read in samfile:
            if read.is_unmapped:
                continue # skip unmap reads
            s = read.query_sequence
            q = ''.join([chr(c+33) for c in read.query_qualities])

            if read.is_reverse:
                s = revComp(s)
                q = q[::-1]

            if ftype == 'fa':
                seq = '>{}\n{}\n'.format(read.query_name, s)
            else:
                seq = '@{}\n{}\n+\n{}\n'.format(read.query_name, s, q)
            of.write(seq)

        of.close()
    
    # output
        

def split_bam_length(bam, output, min=23, max=29):
    """
    split bam file by length
    """
    bam = pysam.AlignmentFile(bam, 'rb')
    fo = pysam.AlignmentFile(output, 'wb', template=bam)

    for read in bam:
        pass


def filt_fq_length(fx, outdir, min=23, max=29, gzipped=True):
    """
    Filter fastq by length
    """
    fname = filename(fx)
    ftype = fx_type(fx)
    fout = '{}.{}'.format(fname, ftype)
    if gzipped:
        fout += '.gz' # gzipped

    # subdir
    sub1 = os.path.join(outdir, '{}_{}'.format(min, max))
    sub2 = os.path.join(outdir, 'not_{}_{}'.format(min, max))
    fx1 = os.path.join(sub1, fout)
    fx2 = os.path.join(sub2, fout)
    check_path(sub1)
    check_path(sub2)

    # run
    if file_exists([fx1, fx2]):
        print('file exists, filt_fx_length() skipped: {}'.format(fname))
    else:
        with xopen(fx, 'rt') as r, xopen(fx1, 'wt') as w1, xopen(fx2, 'wt') as w2:
            for n, s, q, m in readfq(r):
                if ftype == 'fa':
                    seq = '>{}\n{}\n'.format(n, s)
                elif ftype == 'fq':
                    seq = '@{}\n{}\n+\n{}\n'.format(n, s, q)
                else:
                    continue #

                w = w1 if len(s) in range(min, max+1) else w2
                w.write(seq)

    return [fx1, fx2]


def align(fq, index, outdir, k=1, threads=4, gzipped=True, rm_tmp=False):
    """
    Align reads to index
    
    could be: -k1
    """
    fname = filename(fq)
    ftype = fx_type(fq) # fq/fa
    check_path(outdir)

    # files
    bam = os.path.join(outdir, '{}.bam'.format(fname))
    log = os.path.join(outdir, '{}.bowtie.log'.format(fname))
    map_fq = os.path.join(outdir, '{}.{}'.format(fname, ftype))
    unmap_fq = os.path.join(outdir, '{}.unmap.{}'.format(fname, ftype))
    para = '-f' if ftype == 'fa' else '-q'

    # unique + multiple (-k 1)
    cmd = ' '.join([
        'bowtie {} -k {} -S -v 2'.format(para, k),
        '--best -p {}'.format(threads),
        '--un {}'.format(unmap_fq),
        '{} {}'.format(index, fq),
        '2> {}'.format(log),
        '| samtools view -bhS -',
        '| samtools sort -o {} -'.format(bam)])

    cmd_sh = os.path.join(outdir, 'cmd.sh')
    with open(cmd_sh, 'wt') as w:
        w.write(cmd + '\n')

    # run
    if os.path.exists(bam):
        print('file exixts, align() skipped: {}'.format(bam))
    else:
        os.system(cmd)

    # output
    if gzipped:
        map_fq_gz = map_fq + '.gz'
        unmap_fq_gz = unmap_fq + '.gz'
        if file_exists(unmap_fq):
            gzip_cmd(unmap_fq, unmap_fq_gz, decompress=False)
        # convert to fq
        bam_to_fq(bam, map_fq_gz)
        # new
        map_fq_out, unmap_fq_out = [map_fq_gz, unmap_fq_gz]
    else:
        bam_to_fq(bam, map_fq)
        map_fq_out, unmap_fq_out = [map_fq, unmap_fq]

    # output
    return [bam, map_fq_out, unmap_fq_out] # return gzipped files


def align_uniq(fq, index, outdir, threads=4, gzipped=True, rm_tmp=False):
    """
    Align reads to index
    extract unique reads (-m 1)

    unique: -m 1
    multiple: -k 2, ids -> bam
    """
    fname = filename(fq)
    ftype = fx_type(fq) # fq/fa
    check_path(outdir)

    # files
    bam = os.path.join(outdir, '{}.bam'.format(fname))
    log = os.path.join(outdir, '{}.bowtie.log'.format(fname))
    map_fq = os.path.join(outdir, '{}.{}'.format(fname, ftype))
    unmap_fq = os.path.join(outdir, '{}.unmap.{}'.format(fname, ftype))
    para = '-f' if ftype == 'fa' else '-q'

    ##-------------------##
    ## unique reads only
    # unique (-m 1) 
    cmd = ' '.join([
        'bowtie {} -m 1 -S -v 2'.format(para),
        '--un {} --best -p {}'.format(unmap_fq, threads),
        '{} {}'.format(index, fq),
        '2> {}'.format(log),
        '| samtools view -bhS -',
        '| samtools sort -o {} -'.format(bam)])

    # save
    cmd_sh = os.path.join(outdir, 'cmd.sh')
    with open(cmd_sh, 'wt') as w:
        w.write(cmd + '\n')

    # run
    if os.path.exists(bam):
        print('file exists, align skipped: {}'.format(bam))
    else:
        os.system(cmd)

    # output
    if gzipped:
        map_fq_gz = map_fq + '.gz'
        unmap_fq_gz = unmap_fq + '.gz'
        if file_exists(unmap_fq):
            gzip_cmd(unmap_fq, unmap_fq_gz, decompress=False)
        # convert to fq
        bam_to_fq(bam, map_fq_gz)
        # new
        map_fq_out, unmap_fq_out = [map_fq_gz, unmap_fq_gz]
    else:
        bam_to_fq(bam, map_fq)
        map_fq_out, unmap_fq_out = [map_fq, unmap_fq]

    return [bam, map_fq_out, unmap_fq_out] # gzipped output


def align_multi(fq, index, outdir, threads=4, gzipped=True):
    """
    Align reads to index
    extract multiple reads (-k 2)

    multiple: -k 2, ids -> bam
    """
    fname = filename(fq)
    ftype = fx_type(fq) # fq/fa
    check_path(outdir)

    ##-------------------##
    ## unique + multi reads
    both_bam = os.path.join(outdir, '{}.unique_multi_k2.bam'.format(fname))
    both_log = os.path.join(outdir, '{}.unique_multi_k2.bowtie.log'.format(fname))
    para = '-f' if ftype == 'fa' else '-q'

    # cmd
    both_cmd = ' '.join([
        'bowtie {} -k 2 -S -v 2'.format(para),
        '--no-unal --best -p {}'.format(threads),
        '{} {}'.format(index, fq),
        '2> {}'.format(both_log),
        '| samtools view -bhS -',
        '| samtools sort -o {} -'.format(both_bam)])

    # save
    cmd_sh = os.path.join(outdir, 'cmd.sh')
    with open(cmd_sh, 'wt') as w:
        w.write(both_cmd + '\n')

    # run
    if os.path.exists(both_bam):
        print('file exists, align skipped: {}'.format(both_bam))
    else:
        os.system(both_cmd)    

    ##-------------------##
    ## extract multi fq
    multi_ids = os.path.join(outdir, '{}.multi.id.txt'.format(fname))
    multi_fq = os.path.join(outdir, '{}.{}'.format(fname, ftype))
    unmap_fq = os.path.join(outdir, '{}.unmap.{}'.format(fname, ftype))
    multi_fq_gz = multi_fq + '.gz'
    unmap_fq_gz = unmap_fq + '.gz'

    ## cmd
    get_multi_cmd = ' '.join([
        'samtools view -F 0x4 {}'.format(both_bam),
        '| cut -f 1 | sort | uniq -c',
        '| awk \'$1>1 {print $2}\'',
        '> {}'.format(multi_ids),
        '&& seqkit grep -n -f {}'.format(multi_ids),
        '{} > {}'.format(fq, multi_fq),
        '&& seqkit grep -v -n -f {}'.format(multi_ids),
        '{} > {}'.format(fq, unmap_fq)])
    
    # save
    cmd_multi_sh = os.path.join(outdir, 'cmd_multi.sh')
    with open(cmd_multi_sh, 'wt') as w:
        w.write(get_multi_cmd + '\n')

    # run
    if file_exists(multi_fq) or file_exists(multi_fq_gz):
        print('file exists, align skipped: {}'.format(multi_fq))
    else:
        os.system(get_multi_cmd)

    # gzip output
    if gzipped:
        if file_exists(multi_fq): # switch
            gzip_cmd(multi_fq, multi_fq_gz, decompress=False)

        if file_exists(unmap_fq):
            gzip_cmd(unmap_fq, unmap_fq_gz, decompress=False)

        multi_fq_out, unmap_fq_out = [multi_fq_gz, unmap_fq_gz]
    else:
        multi_fq_out, unmap_fq_out = [multi_fq, unmap_fq]

    ##-------------------##
    ## multi alignment/multi_bam
    ## for bam/ -k 1 # only 1 alignment
    multi_bam = os.path.join(outdir, '{}.bam'.format(fname))

    if not file_exists(multi_bam):
        tmp_dir = os.path.join(outdir, 'tmp')
        tmp_bam, _, _ = align(multi_fq_out, index, tmp_dir, k=1, threads=threads, gzipped=gzipped) # 100% mapped
        shutil.copy(tmp_bam, multi_bam)
        shutil.rmtree(tmp_dir)

    return [multi_bam, multi_fq_out, unmap_fq_out]


def split_bam_1u10a(bam, outdir):
    """
    Annotation of input bam file

    1. (sam|bam) extract 1U, 10A reads (1U_10A, 1U_not_10A, not_1U_10A, not_1U_not_10A)     
    """
    fname = filename(bam)
    check_path(outdir)

    # split BAM
    bam_1u_10a = os.path.join(outdir, fname + '.1U_10A.bam')
    bam_1u_10ax = os.path.join(outdir, fname + '.1U_not_10A.bam')
    bam_1ux_10a = os.path.join(outdir, fname + '.not_1U_10A.bam')
    bam_1ux_10ax = os.path.join(outdir, fname + '.not_1U_not_10A.bam')

    # run
    if file_exists([bam_1u_10a, bam_1u_10ax, bam_1ux_10a, bam_1ux_10ax]):
        print('file exists, split_bam_1u10a skipped: {}'.format(fname))
    else:
        if not file_exists(bam + '.bai'):
            pysam.index(bam)
        samfile = pysam.AlignmentFile(bam, 'rb')
        w1 = pysam.AlignmentFile(bam_1u_10a, 'wb', template=samfile)
        w2 = pysam.AlignmentFile(bam_1u_10ax, 'wb', template=samfile)
        w3 = pysam.AlignmentFile(bam_1ux_10a, 'wb', template=samfile)
        w4 = pysam.AlignmentFile(bam_1ux_10ax, 'wb', template=samfile)
        
        for read in samfile.fetch():
            s = read.query_sequence
            q = ''.join([chr(c+33) for c in read.query_qualities])
            if read.is_reverse:
                s = revComp(s)
                q = q[::-1]

            if s[0] == 'T':
                if s[9] == 'A':
                    w1.write(read)
                else:
                    w2.write(read)
            else:
                if s[9] == 'A':
                    w3.write(read)
                else:
                    w4.write(read)

    # output
    return [bam_1u_10a, bam_1u_10ax, bam_1ux_10a, bam_1ux_10ax]
                    

def split_fq_1u10a(fq, outdir=None, gzipped=True):
    """
    Split fastq file by 1U, 10A
    """    
    fname = filename(fq)
    ftype = fx_type(fq)

    # output dir
    if outdir is None:
        outdir = os.path.join(os.path.dirname(fq), '1U_10A')
    check_path(outdir)

    # output
    fq_1u_10a = os.path.join(outdir, '{}.1U_10A.{}'.format(fname, ftype))
    fq_1u_10ax = os.path.join(outdir, '{}.1U_not_10A.{}'.format(fname, ftype))
    fq_1ux_10a = os.path.join(outdir, '{}.not_1U_10A.{}'.format(fname, ftype))
    fq_1ux_10ax = os.path.join(outdir, '{}.not_1U_not_10A.{}'.format(fname, ftype))
    # gzip output
    if gzipped:
        fq_1u_10a += '.gz'
        fq_1u_10ax += '.gz'
        fq_1ux_10a += '.gz'
        fq_1ux_10ax += '.gz'

    fq_list = [fq_1u_10a, fq_1u_10ax, fq_1ux_10a, fq_1ux_10ax]

    if file_exists(fq_list):
        print('file exists, split_fq_1u10a() skipped: {}'.format(fname))
    else:

        fh = xopen(fq, 'rt')
        w1 = xopen(fq_1u_10a, 'wt')
        w2 = xopen(fq_1u_10ax, 'wt')
        w3 = xopen(fq_1ux_10a, 'wt')
        w4 = xopen(fq_1ux_10ax, 'wt')

        for n, s, q, m in readfq(fh):
            if ftype == 'fa':
                seq = '>{}\n{}\n'.format(n, s)
            else:
                seq = '@{}\n{}\n+\n{}\n'.format(n, s, q)

            # check
            if s[0] == 'T':
                if s[9] == 'A':
                    w1.write(seq)
                else:
                    w2.write(seq)
            else:
                if s[9] == 'A':
                    w3.write(seq)
                else:
                    w4.write(seq)

        # close file
        fh.close()
        w1.close()
        w2.close()
        w3.close()
        w4.close()

    return fq_list


def anno_bam(bam, outdir, genome='dm6'):
    """
    Annotate bam alignment, annotationPeaks.pl
    """
    fname = filename(bam)
    check_path(outdir)

    # exe
    anno = shutil.which('annotatePeaks.pl')
    if not anno:
        print('anno_bam() skipped. annotatePeaks.pl not found in $PATH, install HOMER, add to PATH')
        return None

    # check input
    if bam.endswith('.bed'):
        bed = bam # not need
    else:
        # bed = os.path.join(outdir, fname + '.bed')
        bed = os.path.splitext(bam)[0] + '.bed'
        if not file_exists(bed):
            bam_to_bed(bam, bed)
    
    # bed to anno
    anno_stat = os.path.join(outdir, fname + '.anno.stat')
    anno_txt = os.path.join(outdir, fname + '.anno.txt')
    anno_log = os.path.join(outdir, fname + '.anno.log')

    # run
    anno_cmd = ' '.join([
        '{} {} {}'.format(anno, bed, genome),
        '-annStats {}'.format(anno_stat),
        '1> {}'.format(anno_txt),
        '2> {}'.format(anno_log)])

    cmd_sh = os.path.join(outdir, 'cmd_anno.sh')
    with open(cmd_sh, 'wt') as w:
        w.write(anno_cmd + '\n') 

    if file_exists(anno_txt):
        print('file exists, anno_bam() skipped: {}'.format(fname))
    else:
        os.system(anno_cmd)

    # combine anno + seq


def merge_bed_byname(b1, b2, output=None, how='inner'):
    """
    Merge two bed files, by name

    b1 bed file (bed6)
    b2 bed file (bed6)
    output if None, return the pd.DataFrame
    how 'inner', 'outer', 'left', 'right' 
    """
    df1 = pybedtools.BedTool(b1).to_dataframe()
    df2 = pybedtools.BedTool(b2).to_dataframe()
    # merge
    df = pd.merge(df1, df2, left_on='name', right_on='name', how=how)

    if output is None:
        return df
    else:
        df.to_csv(output, index=False, header=False)


def pipe_init(fq, outdir):
    """
    Copy data to dir
    fq_stat()
    """
    fname = filename(fq)
    check_path(outdir)

    fq_to = os.path.join(outdir, os.path.basename(fq))
    if not file_exists(fq_to):
        os.symlink(fq, fq_to)
    fx_stat(fq_to)

    # 1U 10A
    fq_list = split_fq_1u10a(fq_to)
    for f in fq_list:
        fx_stat(f)

    # wrap stat
    df = wrap_stat(outdir)

    return [fq_to, df]


def pipe_collapse(fq, outdir, gzipped=True):
    """
    Collapse, by sequence
    """
    fname = filename(fq)
    check_path(outdir)

    fq_out = fx_collapse(fq, outdir, gzipped=True)
    fx_stat(fq_out)

    # 1U 10A
    fq_list = split_fq_1u10a(fq_out)
    for f in fq_list:
        fx_stat(f)

    # wrap stat
    df = wrap_stat(outdir)

    return [fq_out, df]


def pipe_smRNA(fq, index, outdir, threads=4, genome='dm6', gzipped=True):
    """
    run reads to small RNAs

    return unmap fq
    """
    fname = filename(fq)
    ftype = fx_type(fq)
    check_path(outdir)

    # align
    bam, map_fq, unmap_fq = align(fq, index, outdir, k=1, threads=threads, gzipped=gzipped)

    # annotation
    bam_to_bed(bam)
    anno_bam(bam, outdir, genome)

    # fx stat
    fx_stat(map_fq)

    # 1U 10A
    sub_dir = os.path.join(outdir, '1U_10A')
    bam_list = split_bam_1u10a(bam, sub_dir)
    for b in bam_list:
        f = os.path.splitext(b)[0] + '.' + ftype
        if gzipped:
            f += '.gz'
        bam_to_fq(b, f)
        fx_stat(f)

    # wrap stat
    df = wrap_stat(outdir)

    return [unmap_fq, df]


def pipe_filt_length(fq, outdir, min=23, max=29, gzipped=True):
    """
    filt fastq by length
    """
    fname = filename(fq)
    fq1, fq2 = filt_fq_length(fq, outdir, min=23, max=29, gzipped=gzipped)

    ################
    ## sub-01: 23_29
    ## split fq
    fx_stat(fq1)
    fq1_dir = os.path.dirname(fq1)
    fq1_list = split_fq_1u10a(fq1)
    for i in fq1_list:
        fx_stat(i)

    ## wrap stat
    wrap_stat(fq1_dir)

    #####################
    ## sub-02: not_23_29
    ## split fq
    fx_stat(fq2)
    fq2_dir = os.path.dirname(fq2)
    fq2_list = split_fq_1u10a(fq2)
    for i in fq2_list:
        fx_stat(i)

    ## wrap stat
    wrap_stat(fq2_dir)

    return [fq1, fq2]


def pipe_te(fq, index, outdir, threads=4, genome='dm6', gzipped=True):
    """
    align piRNA reads to TE
    """    
    fname = filename(fq)
    check_path(outdir)

    ##---------------------------------------------------------##
    ## uniq + multiple
    both_dir = os.path.join(outdir, 'unique_multi')
    check_path(both_dir)
    both_bam, both_fq, both_unmap = align(fq, index, both_dir, k=1, threads=threads, gzipped=gzipped)

    # annotation
    bam_to_bed(both_bam)
    anno_bam(both_bam, both_dir, genome)
    both_fq_stat = fx_stat(both_fq)
    both_fq_list = split_fq_1u10a(both_fq, gzipped=gzipped)
    both_bam_list = split_bam_1u10a(both_bam, os.path.join(both_dir, '1U_10A'))
    [fx_stat(f) for f in both_fq_list]

    # wrap read count
    wrap_stat(both_dir)

    ##---------------------------------------------------------##
    ## uniq only
    uniq_dir = os.path.join(outdir, 'unique')
    check_path(uniq_dir)
    uniq_bam, uniq_fq, _ = align_uniq(fq, index, uniq_dir, threads=4, gzipped=gzipped)

    # annotation
    bam_to_bed(uniq_bam)
    anno_bam(uniq_bam, uniq_dir, genome)
    uniq_fq_stat = fx_stat(uniq_fq)
    uniq_fq_list = split_fq_1u10a(uniq_fq, gzipped=gzipped)
    uniq_bam_list = split_bam_1u10a(uniq_bam, os.path.join(uniq_dir, '1U_10A'))
    [fx_stat(f) for f in uniq_fq_list]

    # wrap read count
    wrap_stat(uniq_dir)

    ##---------------------------------------------------------##
    # multi only
    multi_dir = os.path.join(outdir, 'multi')
    check_path(multi_dir)
    multi_bam, multi_fq, _ = align_multi(fq, index, multi_dir, threads=4, gzipped=gzipped)

    # annotation
    bam_to_bed(multi_bam)
    anno_bam(multi_bam, multi_dir, genome)
    multi_fq_stat = fx_stat(multi_fq)
    multi_fq_list = split_fq_1u10a(multi_fq, gzipped=gzipped)
    multi_bam_list = split_bam_1u10a(multi_bam, os.path.join(multi_dir, '1U_10A'))
    [fx_stat(f) for f in multi_fq_list]

    # wrap read count
    wrap_stat(multi_dir)

    return both_unmap


def pipe_overlap(dirA, dirB, outdir):
    """
    Calculate the overlap between TE, piRNA clusters
    subdir: unique, multi, unique_multi
    """
    # b1 = listfile(dirA, "*.bed", recursive=True)
    # b2 = listfile(dirB, "*.bed", recursive=True)

    for sub in ['unique', 'unique_multi', 'multi']:
        # run
        b1 = listfile(os.path.join(dirA, sub), "*.bed")
        b2 = listfile(os.path.join(dirB, sub), "*.bed")

        if len(b1) == 0 or len(b2) == 0:
            log.warning('bed not found, {}'.format(sub))
        else:
            b1 = b1[0] # first one
            b2 = b2[0] # first one

            # output file
            sub_outdir = os.path.join(outdir, sub)
            sub_fn = os.path.basename(b1)
            sub_f = os.path.join(sub_outdir, sub_fn)

            if file_exists(sub_f):
                log.info('file exists, overlap() skipped: {}'.format(sub))
                continue

            # merge two BED file
            df = merge_bed_byname(b1, b2, how='outer')

            # extract count from name
            # id-count
            dm = df['name'].str.split('-', n=1, expand=True)
            if dm.shape[1] == 2:
                # dm = dm.astype({0: str, 1: np.int64})
                dm_chk = dm[1].apply(lambda x:x is None).sum()
                df['count'] = dm[1] if dm_chk == 0 else 1
                df.astype({'count': np.int64})
            else:
                df['count'] = 1

            # re-arrange
            df = df.astype({'count': np.int64}) # numeric
            df.sort_values(by=['count'], ascending=False, inplace=True)
            # df.reset_index(inplace=True)

            # save to file
            check_path(sub_outdir)
            df.to_csv(sub_f, '\t', header=False, index=False)


def read_seqkit_stat(x):
    """
    Read number of reads, in seqkit stat output,
    eg:
    file   format  type  num_seqs    sum_len  min_len  avg_len  max_len
    in.fq  FASTQ   DNA    100,000  2,365,144       18     23.7       36
    """
    df = pd.read_csv(x, delim_whitespace=True, thousands=',')
    df.file = [filename(i) for i in df.file.to_list()]
    # return df.set_index('file').to_dict()
    return df


def get_fq_stat(x):
    """
    return the seqkit stat for fq files
    remove ('.anno.stat')

    u1a10, for 1U_10A files
    """
    l = listfile(x, '*.stat')

    if len(l) > 0:
        l = [read_seqkit_stat(i) for i in l if not 'anno.stat' in i]
        df = pd.concat(l, axis=0)
        df = df[['file', 'num_seqs']]
    else:
        df = pd.DataFrame(columns=['file', 'format', 'num_seqs'])

    # remove 'format' column
    df.drop_duplicates(inplace=True) # fa,fq both exists, 

    return df[['file', 'num_seqs']]


def wrap_stat(outdir):
    """
    wrap reads stat, and 1U_10A
    save in table
    total, unique, multi, ...

    uniq/multi,...
    """
    gname = os.path.basename(outdir.rstrip('/'))
    p = re.compile('^[0-9]{2}\.')
    # if gname in ['unique', 'multi', 'unique_multi']:
    if not p.search(gname):
        t1 = os.path.dirname(outdir.rstrip('/'))
        gname = os.path.basename(t1) + '_' + gname

    # wrap stat
    stat_txt = os.path.join(outdir, 'stat.txt')

    # total
    df1 = get_fq_stat(outdir)
    df1.drop_duplicates(inplace=True) # remove dup
    fname = list(set(df1['file'].to_list()))
    # check
    if len(fname) == 0:
        log.warning('fq.stat not found: {}'.format(outdir))
        return None # break
    else:
        fname = fname[0] # blank
    df1.drop(columns=['file'], inplace=True)
    df1.columns = ['total']
    df1.index = [gname]

    # 1U 10A
    df2 = get_fq_stat(os.path.join(outdir, '1U_10A'))
    df2['file'] = [i.replace(fname + '.', '') for i in df2['file'].to_list()]
    df2x = df2.pivot_table(columns='file', values='num_seqs')
    df2x.index = [gname]

    ## output
    df = pd.concat([df1, df2x], axis=1)
    df.reset_index(inplace=True)

    ## save to file
    df.to_csv(stat_txt, index=False)

    return df


def pipe(fq, outdir, collapse=True, threads=4, genome='dm6', 
    smRNA_index=None, 
    miRNA_index=None,
    te_index=None,
    piRNAcluster_index=None,
    genome_index=None, **kwargs):
    """
    pipeline for piRNA analysis

    length distribution
    alignment stat
    1U, 10A stat
    TE mapping
    piRNA cluster mapping

    overlap TE/piRNA_cluster
    
    Overview:
    1. remove structural RNAs (tRNAs, rRNAs) (unique + multiple)     
    2. filt by length, 23-29 nt   
    3. collapse reads, only compare 1-23 nt (*)  
       collapse reads, (regular)  
       trim reads to 23nt from 3' end
    4. map to TE consensus (unique, multiple, unique + multiple)   
    5. map to piRNA clusters (unique, multiple, unique + multiple)   
    6. map to genome (not-TE, no-piRNAcluster)    
    7. overall, map to genome
    """
    fname = filename(fq)
    check_path(outdir)

    # absolute path
    fq = os.path.abspath(fq)

    # 00 total (symlink)
    log.info('0.copy file: {}'.format(fname))
    dir_00 = os.path.join(outdir, '00.total')
    fq_00, df_00 = pipe_init(fq, dir_00)

    # 01 collapse ()
    log.info('1.collapse')
    if collapse:
        dir_01 = os.path.join(outdir, '01.collapse')
        fq_01, df_01 = pipe_collapse(fq_00, dir_01)
    else:
        log.info('skipped, (if required, use -c, to collapse reads)')
        fq_01 = fq_00

    # 02 to small RNA
    log.info('2.map to small RNA')
    dir_02 = os.path.join(outdir, '02.to_smRNA')
    fq_02, df_02 = pipe_smRNA(fq_01, smRNA_index, dir_02, threads, genome)

    # sys.exit('exit 02')

    # 03 to miRNA
    log.info('3.map to miRNA')
    dir_03 = os.path.join(outdir, '03.to_miRNA')
    fq_03, df_03 = pipe_smRNA(fq_02, miRNA_index, dir_03, threads, genome)

    # 04 filt by length
    log.info('4.filt by length')
    dir_04 = os.path.join(outdir, '04.filt_length')
    fq_04, _ = pipe_filt_length(fq_03, dir_04, min=23, max=29)

    # sys.exit('exit 04')

    # 05 to TE
    log.info('5.map to TE')
    dir_05 = os.path.join(outdir, '05.to_TE')
    fq_05 = pipe_te(fq_04, te_index, dir_05, threads, genome)

    # sys.exit('exit 05')

    # 06 to genome (TE)
    log.info('6.map to genome (TE)')
    dir_06 = os.path.join(outdir, '06.to_genome_TE')
    fq_06 = pipe_te(fq_05, genome_index, dir_06, threads, genome)

    # 07 to piRNA_cluster
    log.info('7.map to piRNA cluster')
    dir_07 = os.path.join(outdir, '07.to_piRNAcluster')
    fq_07 = pipe_te(fq_04, piRNAcluster_index, dir_07, threads, genome)

    # 08 to genome (piRNA cluster)
    log.info('8.map to genome (piRNA cluster)')
    dir_08 = os.path.join(outdir, '08.to_genome_piRNAcluster')
    fq_08 = pipe_te(fq_07, genome_index, dir_08, threads, genome)

    # 09 to genome (all)
    log.info('9.map to genome (all)')
    dir_09 = os.path.join(outdir, '09.to_genome')
    fq_09 = pipe_te(fq_04, genome_index, dir_09, threads, genome)

    # 10 te overlap piRNA cluster
    log.info('10.te_overlap_piRNAcluster')
    dir_10 = os.path.join(outdir, '10.te_overlap_piRNAcluster')
    pipe_overlap(dir_05, dir_07, dir_10)

    # 11 remove temp files
    log.info('11.remove_temp_files')


def main():
    # args = get_args()
    args = vars(get_args())

    ## default
    index_dir = '/home/wangming/data/genome/{}/bowtie_index'.format(args['genome'])
    init_args = {
        'threads': 4,
        'genome': 'dm6',
        'smRNA_index': os.path.join(index_dir, 'smRNA'),
        'miRNA_index': os.path.join(index_dir, 'hairpin'),
        'te_index': os.path.join(index_dir, 'te'),
        'piRNAcluster_index': os.path.join(index_dir, 'piRNA_cluster'),
        'genome_index': os.path.join(index_dir, args['genome'])}

    ## update
    args.update(init_args)

    pipe(fq=args['input'], outdir=args['output'], **args)


if __name__ == '__main__':
    main()

## EOF