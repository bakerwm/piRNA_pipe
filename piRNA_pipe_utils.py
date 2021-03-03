#!/usr/bin/env python3

"""
piRNA analysis (small RNAseq) 


date: 2020-12-23
Author: Ming Wang

# flowchart of piRNA analysis
1. remove structural RNAs (uniqe + multi)
2. remove miRNAs (unique + multi) 
3. remove reads not in [23-29] nt 
4. collapse reads: consider 1-23nt only, allow 1-2 at 3' differ (2019-11-26)
5. split into 4 groups: (1U+/-,10A+/-;)
6. map to TE consensus, (unique + multi), only 1U_not_10A  
7. map to genome (unique + multi) 
Functions:
rename fastq reads: piR0000001-0000001
piR (piRNA), (piRNA number) {reads number}
mapping, unique + multiple/unique



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
from hiseq.utils.helper import * 
from hiseq.utils.seq import Fastx



class OverlapFq(object):
    """Check fastq overlap
    query, 
    subject,
    outdir,
    
    map query to subject using bowtie, check mutations by 'MD' tag
    
    1. convert subject to fasta / collapse
    2. trim to 30 nt
    3. build bowtie-index
    4. alignment
    5. filter bam MD: tag in sam file
    6. return the query sequence
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'query': None,
            'subject': None,
            'outdir': None,
            'mm': 2,
            'mm_range': [5, 28],
            'threads': 4
        }
        self = update_obj(self, args_init, force=False)

        if self.query == None or self.subject == None:
            raise Exception('query, subject required')

        if self.outdir == None:
            self.outdir = str(pathlib.Path.cwd())

        # files
        s_dir = os.path.join(self.outdir, 'index')
        if self.is_bowtie_index(self.subject):
            self.s_name = os.path.basename(self.subject)
            self.s_index = self.subject
        else:
            self.s_name = fq_name(self.subject)
            self.s_fa = os.path.join(s_dir, self.s_name+'.fa')
            self.s_index = os.path.join(s_dir, self.s_name)
        # other files
        self.q_name = fq_name(self.query)
        self.q_bam = os.path.join(self.outdir, self.q_name+'.align.bam')
        self.q_bam_filt = self.q_bam.replace('.bam', '.filt.bam')
        self.q_align_log = os.path.join(self.outdir, self.q_name+'.align.log')
        self.out_fq_plain = os.path.join(self.outdir, self.q_name+'.fq') # gzip
        self.out_fq = self.out_fq_plain + '.gz'

        # dirs
        check_path(s_dir)


    def is_bowtie_index(self, x):
        if isinstance(x, str):
            return file_exists(x+'.1.ebwt')


    def prep_index(self):
        """Collapse, index
        bowtie-build
        """
        if self.is_bowtie_index(self.subject):
            log.info('prep_index() skipped, subject is bowtie index')
        else:
            #1.collapse
            if not file_exists(self.s_fa):
                Fastx(self.subject).collapse(self.s_fa, fq_out=False)

            #2.build bowtie index
            cmd = ' '.join([
                'bowtie-build -q',
                '--threads {}'.format(self.threads),
                '{}'.format(self.s_fa),
                '{}'.format(self.s_index)
            ])
            if not self.is_bowtie_index(self.s_index):
                run_shell_cmd(cmd)


    def run_align(self):
        """Align query to subject
        bowtie -v 2 -S -k 1 -v 2 -l 23 -n 3
        """
        cmd = ' '.join([
            'bowtie -S -k 1 -v 2',
            '-p {}'.format(self.threads),
            '-l 23 -n 3',
            '--no-unal',
            '-x {}'.format(self.s_index),
            '{}'.format(self.query),
            '2>{}'.format(self.q_align_log),
            '| samtools view -bhS -',
            '| samtools sort -o {} -'.format(self.q_bam),
            '&& samtools index {}'.format(self.q_bam)
        ])
        if not file_exists(self.q_bam):
            run_shell_cmd(cmd)


#     def check_align(self, x):
#         """Check MD in sam record
#         make sure MD:Z:20G0G6 
#         mutations at 5-28 ranges
#         mm <=3
#         x AlignedSegment
#         """
#         f = False
#         if isinstance(x, pysam.AlignedSegment):
#             if x.qstart == 0 and x.qstart == x.reference_start:
#                 md = [i[1] for i in x.tags if i[0] == 'MD'] # ['28G0G1'], ['29']
#                 if len(md) == 1:
#                     # mismatches ?!
#                     mm = re.findall('[ACGT]', md[0])
#                     p = re.split('[ACGTN]', md[0]) # '[^0-9]'
#                     if len(mm) <= self.mm:
#                         f = True # mismatches
#                     if len(p) > 1:
#                         px = sum(list(map(int, p))[:-1]) # remove last one
#                         f = px in range(self.mm_range[0], self.mm_range[1])
#                     else:
#                         f = True # ?
#         return f


    def check_align(self, x):
        """Check MD in sam record
        make sure MD:Z:20G0G6 
        mutations at 5-28 ranges
        mm <=3
        x AlignedSegment
        """
        f = False
        if isinstance(x, pysam.AlignedSegment):
            print(x.qstart, x.reference_start)
            if x.qstart == 0 and x.qstart == x.reference_start:
                md = [i[1] for i in x.tags if i[0] == 'MD'] # ['28G0G1'], ['29']
                if len(md) == 1:
                    # mismatches ?!
                    print('AAAA-1')
                    mm = re.findall('[ACGT]', md[0])
                    p = re.split('[ACGTN]', md[0]) # '[^0-9]'
                    if len(mm) <= self.mm:
                        print('AAAA-2')
                        f = True # mismatches
                    if len(p) > 1:
                        print('AAAA-3')
                        px = sum(list(map(int, p))[:-1]) # remove last one
                        f = px in range(self.mm_range[0], self.mm_range[1])
                    else:
                        f = True # ?
        return f


    def run_filter(self):
        """Filt alignments, mismatches within range[]
        pysam.AlignmentFiles()
        """
        samfile = pysam.AlignmentFile(self.q_bam, 'rb')
        if not file_exists(self.q_bam_filt):
            destfile = pysam.AlignmentFile(self.q_bam_filt, 'wb', 
                template=samfile)
            for r in samfile:
                if self.check_align(r):
                    destfile.write(r)


    def save_to_fq(self):
        """Save filt bam file, to fastq format
        out_fq
        """
        if not file_exists(self.out_fq_plain):
            cmd = ' '.join([
                'samtools fastq',
                '-@ {}'.format(self.threads),
                '{}'.format(self.q_bam_filt),
                '| gzip > {}'.format(self.out_fq)
            ])
            run_shell_cmd(cmd)
            # pysam.fastq('-@', '{}'.format(self.threads),
            #     '-s', self.out_fq_plain, 
            #     self.q_bam_filt, save_stdout=self.out_fq_plain)

        # if not file_exists(self.out_fq):
        #     gzip_cmd(self.out_fq_plain, self.out_fq, decompress=False)
        
    
    def stat(self):
        """Stat, overlap reads, input reads"""
        query_c = Fastx(self.query).number_of_seq()
        ov_c = Fastx(self.out_fq).number_of_seq()
        ov_pct = ov_c/query_c*100
        # message
        msg = 'query: {}, overlap: {}, percent: {:.2f}%'.format(
            query_c, ov_c, ov_pct
        )
        print(msg)
        return (query_c, ov_c, ov_pct)        


    def run(self):
        self.prep_index()
        self.run_align()
        self.run_filter()
        self.save_to_fq()
        self.stat()
        return self.out_fq

















################################################################################
def collapse_fx(fx_in, fx_out):
    """
    collapse fastx by seq
    fx_out: could be : fastq/fasta, gz
    """
    outdir = os.path.dirname(fx_out)
    check_path(outdir)
    # check fa/fq
    fx_name = os.path.basename(fx_out)
    fx_name = fx_name.replace('.gz', '')
    if fx_name.endswith('.fq') or fx_name.endswith('.fastq'):
        fq_out = True
    else:
        fq_out = False
    if file_exists(fx_out):
        log.info('collapse_fx() skipped, file exists: {}'.format(fx_out))
    else:
        d = {}
        with xopen(fx_in, 'rt') as r:
            for n, s, q, m in readfq(r):
                d[s] = d.get(s, 0) + 1

        # sort by value (decreasing)
        i = 0
        with xopen(fx_out, 'wt') as w:
            for k, v in sorted(d.items(), key=lambda item: item[1], reverse=True):
                i += 1
                out = '>{}-{}\n{}\n'.format(i, v, k)
                if fq_out:
                    out = '@{}-{}\n{}\n+\n{}\n'.format(i, v, k, 'I'*len(k))
                w.write(out)
    return fx_out


def count_fx(fx, collapsed=False):
    """Count reads
    if collapsed
    parse the name: id-count
    return 
    counts:
    species:
    """
    if not isinstance(fx, str):
        log.warning('count_fx() skipped, fx expect str, got {}'.format(
            type(fx).__name__))
    else:
        if not file_exists(fx):
            log.warning('count_fx() skipped, fx not exists: {}'.format(fx))
        else:
            fname = filename(fx)
            out_toml = os.path.join(os.path.dirname(fx), fname+'.fx_stat.toml')
            if file_exists(out_toml):
                log.info('count_fx() skipped, file exists: {}'.format(out_toml))
                d = Toml().from_toml(out_toml)
            else:
                n_read = 0
                n_species = 0
                if collapsed:
                    with xopen(fx, 'rt') as r:
                        for n, _, _, _ in readfq(r):
                            if '-' in n:
                                name, i = n.split('-', 1)
                                i = eval(i.lstrip('0'))
                            else:
                                i = 1
                            n_read += i
                            n_species += 1
                else:
                    n_read = n_species = Fastx(fx).number_of_seq()
                # save to file
                d = {
                    'name': fname, 
                    'num_seqs': n_read, 
                    'num_species': n_species
                }
                Toml(d).to_toml(out_toml)
            # return pd.DataFrame.from_dict(d, 'columns')
            return {fname: d}


def count_fx_dir(x, collapsed=False):
    """Count reads in directory: x
    Concat all the counts
    assume all fx files are gzipped
    """
    d = {}
    if isinstance(x, str):
        if os.path.isdir(x):
            fx_list = listfile(x, '*.gz', recursive=False)
            if len(fx_list) > 0:
                for fx in fx_list:
                    d.update(count_fx(fx, collapsed))
                # save to file
                out_toml = os.path.join(x, 'fx_stat.toml')
                Toml(d).to_toml(out_toml)
        else:
            log.warning('count_fx_dir() skipped, x not a directory: {}'.format(x))


def stat_fq(fq, outdir=None):
    """
    seqkit stat fq
    output: fq.fq_stat.txt
    """
    if isinstance(fq, str):
        if not isinstance(outdir, str):
            outdir = os.path.dirname(fq)
        fname = fq_name(fq)
        fq_stat_name = '{}.fq_stat.txt'.format(fq_name)
        fq_stat_txt = os.path.join(outdir, fq_stat_name)
        if outdir is None:
            outdir = os.path.dirname(fq)
        check_path(outdir)
        # run    
        cmd = 'seqkit stat {} > {}'.format(fq, fq_stat_txt)
        if file_exists(fq_stat_txt):
            # log.info('stat_fq() skipped, file exists: {}'.format(fname))
            pass
        else:
            # log.info('Run stat_fq(): {}'.format(fname))
            os.system(cmd)
    elif isinstance(fq, list):
        if len(fq) > 0:
            [stat_fq(i) for i in fq]
    else:
        log.warning('stat_fq() skipped, str, list expected, got {}'.format(
            type(fq).__name__))


def stat_fq2(x, recursive=False):
    """Run seqkit stat for fastq files in x
    recursive
    """
    if isinstance(x, str):
        if os.path.isdir(x):
            fq_list = list_fq(x, recursive)
            stat_fq(fq_list)
            # organize
            wrap_stat_fq(x, recursive)


def read_seqkit_stat(x):
    """
    Read number of reads, in seqkit stat output,
    eg:
    file   format  type  num_seqs    sum_len  min_len  avg_len  max_len
    in.fq  FASTQ   DNA    100,000  2,365,144       18     23.7       36
    """
    if isinstance(x, str):
        if check_file(x, emptycheck=True):
            df = pd.read_csv(x, delim_whitespace=True, thousands=',')
            df.file = [fq_name(i) for i in df.file.to_list()]
        else: 
            log.warning('read_seqkit_stat() skipped, empty file: {}'.format(x))
            df = None
    else:
        log.warning('read_seqkit_stat() expect single file, skipped ...')
        df = None

    return df


def parse_stat_fq(x):
    """Return the seqkit stat files in directory: x
    all fq_stat.txt files in the same dir of the x file
    """
    s_list = [] # empty
    if isinstance(x, str):
        if os.path.isdir(x):
            s_list = listfile(x, '*.fq_stat.txt')

    # reading files
    if len(s_list) > 0:
        s_frames = list(map(read_seqkit_stat, s_list))
        s_frames = [i for i in s_frames if i is not None]
        if len(s_frames) > 0:
            df = pd.concat(s_frames, axis=0)
            df = df[['file', 'num_seqs']] # select columns
        else:
            log.warning('parse_stat_fq() skipped, empty fq_stat.txt files')
            df = pd.DataFrame(columns=['file', 'num_seqs'])
    else:
        log.warning('parse_stat_fq() skipped, no *.fq_stat.txt found')
        df = pd.DataFrame(columns=['file', 'num_seqs'])

    # remove 'format' column
    df.drop_duplicates(inplace=True) # fa,fq both exists, 

    return df


def wrap_stat_fq(x, recursive=False):
    """wrap reads stat, in directory: x
    level-1: 00.raw_data/*.fq.gz
    """
    if isinstance(x, str):
        if os.path.isdir(x):
            x_dirs = [x]
            # subdirs
            if recursive:
                x_subdirs = [i for i in listdir(x, recursive=True, 
                    include_dir = True) if os.path.isdir(i)]
                x_dirs.extend(x_subdirs)

            # wrap stat
            for d in x_dirs:
                # stat_json = os.path.join(d, 'fq_stat.json')
                stat_toml = os.path.join(d, 'fq_stat.toml')
                df = parse_stat_fq(d)
                df = df.set_index('file').to_dict() # .to_json(stat_json)
                if len(df.keys()) > 0:
                    Toml(df).to_toml(stat_toml)


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


def split_fq_1u10a(fq, outdir=None, gzipped=True, remove=False):
    """Split fastq file by 1U, 10A
    1U   +/-
    10A  +/-
    """
    fname = fq_name(fq)
    ftype = Fastx(fq).format # fasta/fastq
    # output dir
    if outdir is None:
        outdir = os.path.join(os.path.dirname(fq), '1U_10A')
    check_path(outdir)

    # output
    f_ext = ftype + '.gz' if gzipped else ftype # fastq/fasta
    fq_names = ['1U_10A', '1U_not_10A', 'not_1U_10A', 'not_1U_not_10A']
    fq_files = [os.path.join(outdir, '{}.{}.{}'.format(fname, i, f_ext)) \
        for i in fq_names]

    if not check_file(fq, emptycheck=True):
        log.error('split_fq_1u10a() skipped, file is empty: {}'.format(fname))
    elif all(file_exists(fq_files)):
        log.info('split_fq_1u10a() skipped, file exists: {}'.format(fname))
    else:
        fh = xopen(fq, 'rt')
        w1 = xopen(fq_files[0], 'wt')
        w2 = xopen(fq_files[1], 'wt')
        w3 = xopen(fq_files[2], 'wt')
        w4 = xopen(fq_files[3], 'wt')

        for n, s, q, m in readfq(fh):
            if ftype == 'fasta':
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

    # remove 1U 10A files
    if remove:
        file_remove(fq_files, ask=False)

    return fq_files


def split_bam_1u10a(bam, outdir=None):
    """Check 1U10A for reads in bam file

    1. (sam|bam) extract 1U, 10A reads (1U_10A, 1U_not_10A, not_1U_10A, not_1U_not_10A)     
    """
    fname = filename(bam)

    # output dir
    if outdir is None:
        outdir = os.path.join(os.path.dirname(bam), '1U_10A')
    check_path(outdir)

    if not file_exists(bam + '.bai'):
        pysam.index(bam)

    # split BAM
    bam_names = ['1U_10A', '1U_not_10A', 'not_1U_10A', 'not_1U_not_10A']
    bam_files = [os.path.join(outdir, fname + i + '.bam') for i in bam_names]

    if all(file_exists(bam_files)):
        log.info('file exists, split_bam_1u10a skipped: {}'.format(fname))
    else:
        samfile = pysam.AlignmentFile(bam, 'rb')
        w1 = pysam.AlignmentFile(bam_files[0], 'wb', template=samfile)
        w2 = pysam.AlignmentFile(bam_files[1], 'wb', template=samfile)
        w3 = pysam.AlignmentFile(bam_files[2], 'wb', template=samfile)
        w4 = pysam.AlignmentFile(bam_files[3], 'wb', template=samfile)
        
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
    return bam_files
 

def splilt_fq_size(fq, outdir=None, min=23, max=29, gzipped=True):
    """Split fastq file by the size 
    23-29 nt
    """
    fname = filename(fq)
    ftype = Fastx(fq).format

    # output dir
    if outdir is None:
        outdir = os.path.join(os.path.dirname(fq), '1U_10A')
    # outfiles
    f_ext = ftype + '.gz' if gzipped else ftype
    f_name = fname + '.' + f_ext
    dir_in = os.path.join(outdir, '{}_{}'.format(min, max))
    dir_ex = os.path.join(outdir, 'not_{}_{}'.format(min, max))
    f_in = os.path.join(dir_in, f_name)
    f_ex = os.path.join(dir_ex, f_name)
    check_path([dir_in, dir_ex])

    # run
    if not check_file(fq, emptycheck=True):
        log.error('file is empty, split_fq_size() skipped: {}'.format(fname))
    if all(file_exists([f_in, f_ex])):
        log.info('file exists, split_fq_length() skipped: {}'.format(fname))
    else:
        with xopen(fq, 'rt') as r, \
            xopen(f_in, 'wt') as w1, \
            xopen(f_ex, 'wt') as w2:
            for n, s, q, m in readfq(r):
                if ftype == 'fasta':
                    seq = '>{}\n{}\n'.format(n, s)
                elif ftype == 'fastq':
                    seq = '@{}\n{}\n+\n{}\n'.format(n, s, q)
                else:
                    continue #

                w = w1 if len(s) in range(min, max+1) else w2
                w.write(seq)

    return [f_in, f_ex]


def split_bam_size(bam, output=None, min=23, max=29):
    """
    split bam file by length
    """
    fname = filename(fq)    
    ftype = Fastx(fq).format

    # output dir
    if outdir is None:
        outdir = os.path.join(os.path.dirname(fq), '1U_10A')

    bam = pysam.AlignmentFile(bam, 'rb')
    fo = pysam.AlignmentFile(output, 'wb', template=bam)

    for read in bam:
        pass


################################################################################
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
        log.info('file exists, bam_to_bed() skipped: {}'.format(os.path.basename(bam)))
    else:
        # pybedtools.BedTool(bam).bam_to_bed().saveas(bed)
        cmd = 'bedtools bamtobed -i {} > {}'.format(bam, bed)
        os.system(cmd)


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
        fq = os.path.join(os.path.dirname(bam), fname + '.fastq.gz') # default
    # output fmt
    fq_suffix = os.path.basename(fq)
    if fq_suffix.endswith('.gz'):
        fq_suffix = fq_suffix.rstrip('.gz')
    fq_ext = os.path.splitext(fq_suffix)[1]
    fq_format = 'fasta' if fq_ext.lower() in ['.fa', '.fasta'] else 'fastq'

    # read bam
    samfile = pysam.AlignmentFile(bam, 'rb')

    if file_exists(fq):
        log.info('file exists, bam_to_fq() skipped: {}'.format(os.path.basename(fq)))
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

            if fq_format == 'fasta':
                seq = '>{}\n{}\n'.format(read.query_name, s)
            else:
                seq = '@{}\n{}\n+\n{}\n'.format(read.query_name, s, q)
            of.write(seq)

        of.close()


# unique + multi: k=1
def align_both(fq, index, outdir, k=1, threads=4, gzipped=True, rm_tmp=False):
    """
    Align reads to index
    could be: -k1
    """
    fname = filename(fq)
    ftype = Fastx(fq).format # fasta/fastq
    check_path(outdir)
    # files
    bam = os.path.join(outdir, '{}.bam'.format(fname))
    align_log = os.path.join(outdir, '{}.bowtie.log'.format(fname))
    map_fq = os.path.join(outdir, '{}.{}'.format(fname, ftype))
    unmap_fq = os.path.join(outdir, '{}.unmap.{}'.format(fname, ftype))
    para = '-f' if ftype == 'fasta' else '-q'
    # unique + multiple (-k 1)
    cmd = ' '.join([
        'bowtie {} -k {} -S -v 2'.format(para, k),
        '--best -p {}'.format(threads),
        '--un {}'.format(unmap_fq),
        '-x {} {}'.format(index, fq),
        '2> {}'.format(align_log),
        '| samtools view -bhS -',
        '| samtools sort -o {} -'.format(bam),
        '&& samtools index {}'.format(bam)])

    cmd_sh = os.path.join(outdir, 'cmd.sh')
    with open(cmd_sh, 'wt') as w:
        w.write(cmd + '\n')
    # run
    if os.path.exists(bam):
        log.info('file exixts, align() skipped: {}'.format(bam))
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
        map_fq_out, unmap_fq_out = [map_fq_gz, unmap_fq_gz]
    else:
        bam_to_fq(bam, map_fq)
        map_fq_out, unmap_fq_out = [map_fq, unmap_fq]
    # output
    return [bam, map_fq_out, unmap_fq_out] # return gzipped files


# unique: m=1
def align_uniq(fq, index, outdir, threads=4, gzipped=True, rm_tmp=False):
    """
    Align reads to index
    extract unique reads (-m 1)

    unique: -m 1
    multiple: -k 2, ids -> bam
    """
    fname = filename(fq)
    ftype = Fastx(fq).format
    check_path(outdir)
    # files
    bam = os.path.join(outdir, '{}.bam'.format(fname))
    align_log = os.path.join(outdir, '{}.bowtie.log'.format(fname))
    map_fq = os.path.join(outdir, '{}.{}'.format(fname, ftype))
    unmap_fq = os.path.join(outdir, '{}.unmap.{}'.format(fname, ftype))
    para = '-f' if ftype == 'fasta' else '-q'
    ##-------------------##
    # unique (-m 1) 
    cmd = ' '.join([
        'bowtie {} -m 1 -S -v 2'.format(para),
        '--un {} --best -p {}'.format(unmap_fq, threads),
        '-x {} {}'.format(index, fq),
        '2> {}'.format(align_log),
        '| samtools view -bhS -',
        '| samtools sort -o {} -'.format(bam),
        '&& samtools index {}'.format(bam)])
    # save
    cmd_sh = os.path.join(outdir, 'cmd.sh')
    with open(cmd_sh, 'wt') as w:
        w.write(cmd + '\n')
    # run
    if os.path.exists(bam):
        log.info('file exists, align skipped: {}'.format(bam))
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
        map_fq_out, unmap_fq_out = [map_fq_gz, unmap_fq_gz]
    else:
        bam_to_fq(bam, map_fq)
        map_fq_out, unmap_fq_out = [map_fq, unmap_fq]
    return [bam, map_fq_out, unmap_fq_out] # gzipped output


# multi: 
def align_multi(fq, index, outdir, threads=4, gzipped=True):
    """
    Align reads to index
    extract multiple reads (-k 2)

    multiple: -k 2, ids -> bam
    """
    fname = filename(fq)
    ftype = Fastx(fq).format
    check_path(outdir)
    ##-------------------##
    ## unique + multi reads
    both_bam = os.path.join(outdir, '{}.unique_multi_k2.bam'.format(fname))
    both_log = os.path.join(outdir, '{}.unique_multi_k2.bowtie.log'.format(fname))
    para = '-f' if ftype == 'fasta' else '-q'
    # cmd
    both_cmd = ' '.join([
        'bowtie {} -k 2 -S -v 2'.format(para),
        '--no-unal --best -p {}'.format(threads),
        '-x {} {}'.format(index, fq),
        '2> {}'.format(both_log),
        '| samtools view -bhS -',
        '| samtools sort -o {} -'.format(both_bam),
        '&& samtools index {}'.format(both_bam)])
    # save
    cmd_sh = os.path.join(outdir, 'cmd.sh')
    with open(cmd_sh, 'wt') as w:
        w.write(both_cmd + '\n')
    # run
    if os.path.exists(both_bam):
        log.info('file exists, align skipped: {}'.format(both_bam))
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
        log.info('file exists, align skipped: {}'.format(multi_fq))
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
        tmp_bam, _, _ = align_both(multi_fq_out, index, tmp_dir, k=1, 
            threads=threads, gzipped=gzipped) # 100% mapped
        shutil.copy(tmp_bam, multi_bam)
        shutil.rmtree(tmp_dir)
    return [multi_bam, multi_fq_out, unmap_fq_out]


def anno_bam(bam, outdir, genome='dm6'):
    """
    Annotate bam alignment, annotationPeaks.pl
    """
    fname = filename(bam)
    check_path(outdir)

    # exe
    anno = shutil.which('annotatePeaks.pl')
    if not anno:
        log.info('anno_bam() skipped. annotatePeaks.pl not found in $PATH, \
            install HOMER, add to PATH')
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
        log.info('file exists, anno_bam() skipped: {}'.format(fname))
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


def pipe_align(fq, index, outdir, threads=4, genome='dm6', gzipped=True,
    unique_multi='both', remove_1u10a=False):
    """Align reads to unique, multiple 
    unique_multi str unique, multi, both, all
    """
    f_name = fq_name(fq)
    f_name = f_name.replace('.unmap', '') # remove suffix
    fq_align = os.path.join(outdir, f_name + '.fastq.gz')

    ## priority: both > unique > multi
    ##---------------------------------------------------------##
    # multi only
    if unique_multi in ['multi', 'all']:
        multi_dir = os.path.join(outdir, 'multi')
        check_path(multi_dir)
        multi_bam, multi_fq, multi_unmap = align_multi(fq, index, multi_dir, 
            threads=4, gzipped=gzipped)
        # annotation
        anno_bam(multi_bam, multi_dir, genome)
        multi_bam_list = split_bam_1u10a(multi_bam)
        # output
        fq_unmap = multi_unmap
        file_symlink(multi_fq, fq_align)
    
    ##---------------------------------------------------------##
    ## uniq only
    if unique_multi in ['unique', 'all']:
        uniq_dir = os.path.join(outdir, 'unique')
        check_path(uniq_dir)
        uniq_bam, uniq_fq, uniq_unmap = align_uniq(fq, index, uniq_dir, 
            threads=4, gzipped=gzipped)
        # annotation
        anno_bam(uniq_bam, uniq_dir, genome)
        uniq_bam_list = split_bam_1u10a(uniq_bam)
        # output
        fq_unmap = uniq_unmap
        if file_exists(fq_align):
            file_remove(fq_align, ask=False)
        file_symlink(uniq_fq, fq_align)

    ##---------------------------------------------------------##
    ## uniq + multiple
    if unique_multi in ['both', 'all']:
        both_dir = os.path.join(outdir, 'unique_multi')
        check_path(both_dir)
        both_bam, both_fq, both_unmap = align_both(fq, index, both_dir, k=1, 
            threads=threads, gzipped=gzipped)
        # annotation
        anno_bam(both_bam, both_dir, genome)
        both_bam_list = split_bam_1u10a(both_bam)
        # output
        fq_unmap = both_unmap
        if file_exists(fq_align):
            file_remove(fq_align, ask=False)
        file_symlink(both_fq, fq_align)

    return (fq_align, fq_unmap)


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



################################################################################
# Deprecated
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
    stat_fq(fq_to)

    # 1U 10A
    fq_list = split_fq_1u10a(fq_to)
    for f in fq_list:
        stat_fq(f)

    # wrap stat
    df = wrap_stat_fq(outdir)

    return [fq_to, df]


def pipe_collapse(fq, outdir, gzipped=True):
    """
    Collapse, by sequence
    """
    fname = filename(fq)
    check_path(outdir)

    fq_out = collapse_fx(fq, outdir, gzipped=True)
    stat_fq(fq_out)

    # 1U 10A
    fq_list = split_fq_1u10a(fq_out)
    for f in fq_list:
        stat_fq(f)

    # wrap stat
    df = wrap_stat_fq(outdir)

    return [fq_out, df]


def pipe_smRNA(fq, index, outdir, threads=4, genome='dm6', gzipped=True):
    """
    run reads to small RNAs

    return unmap fq
    """
    fname = filename(fq)
    ftype = Fastx(fq).format
    check_path(outdir)

    # align
    bam, map_fq, unmap_fq = align(fq, index, outdir, k=1, threads=threads, gzipped=gzipped)

    # annotation
    bam_to_bed(bam)
    anno_bam(bam, outdir, genome)

    # fq stat
    stat_fq(map_fq)

    # 1U 10A
    sub_dir = os.path.join(outdir, '1U_10A')
    bam_list = split_bam_1u10a(bam, sub_dir)
    for b in bam_list:
        f = os.path.splitext(b)[0] + '.' + ftype
        if gzipped:
            f += '.gz'
        bam_to_fq(b, f)
        stat_fq(f)

    # wrap stat
    df = wrap_stat_fq(outdir)

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
    stat_fq(fq1)
    fq1_dir = os.path.dirname(fq1)
    fq1_list = split_fq_1u10a(fq1)
    for i in fq1_list:
        stat_fq(i)

    ## wrap stat
    wrap_stat_fq(fq1_dir)

    #####################
    ## sub-02: not_23_29
    ## split fq
    stat_fq(fq2)
    fq2_dir = os.path.dirname(fq2)
    fq2_list = split_fq_1u10a(fq2)
    for i in fq2_list:
        stat_fq(i)

    ## wrap stat
    wrap_stat_fq(fq2_dir)

    return [fq1, fq2]




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

def fa2fq(fa, fq):
    """
    Convert fasta to fastq
    "I" for quality
    """
    fq_dir = os.path.dirname(fq)
    check_path(fq_dir)

    if file_exists(fq):
        log.info('file exists, fa2fq() skipped:')
    else:
        with xopen(fa, 'rt') as r, xopen(fq, 'wt') as w:
            for n, s, q, m in readfq(r):
                fq = '@{}\n{}\n+\n{}\n'.format(n, s, q)
                w.write(fq)