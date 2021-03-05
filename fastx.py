#!/usr/bin/env python3

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
import pyfastx
import pysam
import pandas as pd
import numpy as np
from xopen import xopen
from multiprocessing import Pool
from hiseq.utils.helper import * 
from hiseq.utils.seq import Fastx
from utils import get_fx_name



class FxFragSize(object):
    """Calculate length distribution for fx
    fx list of files (could be fastq, bam)
    outdir
    
    Example:
    >>> FxFragSize(fx).run()
    """
    def __init__(self, fx, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.fx = fx
        self.init_args()
        
        
    def init_args(self):
        args_init = {
            'outdir': None,
            'parallel_jobs': 1
        }
        self = update_obj(self, args_init, force=False)
        if isinstance(self.fx, str):
            self.fx = [self.fx]
        if not isinstance(self.fx, list):
            raise ValueError('fx, expect list, got {}'.format(
                type(self.fx).__name__))
        self.fx = [i for i in self.fx if check_file(i, emptycheck=True)]
        if len(self.fx) == 0:
            raise ValueError('fx, not found')


    def run_single(self, i):
        """Stat read size for BAM/fx
        save as .csv file

        >fragsize.csv
        length strand count
        """
        fx = self.fx[i]
        f_name, f_ext = os.path.splitext(os.path.basename(fx))
        # output
        outdir = self.outdir if isinstance(self.outdir, str) else os.path.dirname(fx)
        if f_ext == '.gz':
            f_name = os.path.splitext(f_name)[0]
        csv_file = os.path.join(outdir, f_name + '.fragsize.csv')
        if file_exists(csv_file):
            log.info('fragsize() skipped, file exists: {}'.format(csv_file))
        else:
            if f_ext == '.bam':
                BamFragSize(fx, asPE=False, strandness=True).saveas(csv_file)
            elif f_ext == '.gz':
                Fastx(fx).len_dist(csv_file=csv_file)
            else:
                pass


    def run(self):
        """List all bam files
        *.bam files
        *.fq.gz files
        """
        if len(self.fx) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single, 
                    range(len(self.fx)))
        else:
            for x in range(len(self.fx)):
                self.run_single(x)



################################################################################
class FxCount(object):
    """Count number of fastx records
    parameters:
    fx
    
    output:
    name, seq, unique
        
    for fasta (collapsed), >1-200,
    add the total amount
    
    Example:
    >>> FxCount(fx).run()
    """
    def __init__(self, fx, **kwargs):
        self.fx = fx
        self.init_args()
        
        
    def init_args(self):
        args_init = {
            'outdir': None,
            'parallel_jobs': 4,
        }
        self = update_obj(self, args_init, force=False)
        if isinstance(self.fx, str):
            self.fx = [self.fx]
        # exists
        self.fx = [i for i in self.fx if check_file(i, emptycheck=True)]
        if len(self.fx) == 0:
            raise ValueError('No fastx files')
        
        
    def count_fq(self, x):
        """Return the nmber of records
        seq
        """
#         n_seq = Fastx(x).number_of_seq()
#         return (n_seq, n_seq)  # unique, total    
        f = pyfastx.Fastx(x)
        i = 0
        s = []
        for name,_,_,_ in f:
            i += 1
            try:
                a,b = name.split('-', 1) # make sure, fastx_collapser output
                s.append(int(b.lstrip('0')))
            except:
                s.append(1)
        return (i, sum(s)) # unique, total
        
    
    def count_fa(self, x):
        """Count reads for fasta, with name: 
        
        example:
        >1-200
        AACCAACC
        >2-188
        CCGGTTAAC
        add the total numbers
        """        
        f = pyfastx.Fastx(x)
        i = 0
        s = []
        for name,_,_ in f:
            i += 1
            try:
                a,b = name.split('-', 1) # make sure, fastx_collapser output
                s.append(int(b.lstrip('0')))
            except:
                s.append(1)
        return (i, sum(s)) # unique, total
    
    
    def count(self, x):
        """Count reads for fastx"""
        # x_name = fq_name(x)
        x_name = get_fx_name(x, fix_unmap=False)
        x_type = Fastx(x).format
        # save to file
        x_stat = os.path.join(os.path.dirname(x), x_name+'.fx_stat.toml')
        if file_exists(x_stat):
            log.info('FxCount() skipped, file exists: {}'.format(x))
            stat = Config().load(x_stat)
        else:
            n_seq = self.count_fq(x) if x_type == 'fastq' else self.count_fa(x)
            stat = {
                'name': x_name,
                'type': x_type,
                'unique': n_seq[0],
                'total': n_seq[1]            
            }
            Config().dump(stat, x_stat)
        return stat
        
        
    def run(self):
        """For n files"""
        s = [self.count(i) for i in self.fx]
        # msg
        s_list = ['\t'.join(['fx_name', 'fx_type', 'unique', 'total'])]
        s_list.extend(['\t'.join(map(str, i.values())) for i in s])
        msg = '\n'.join(s_list)
        print(msg)



################################################################################
class FxU1A10(object):
    """Check U1 and A10 for fastx file
    check for fx
    
    Example:
    >>> FxU1A10(fx).run()
    """
    def __init__(self, fx, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.fx = fx
        self.init_args()
       
    
    def init_args(self):
        """Require: fx
        """
        args_init = {
            'outdir': None,
            'parallel_jobs': 4,            
        }
        self = update_obj(self, args_init, force=False)
        # file exists, not empty
        if isinstance(self.fx, str):
            self.fx = [self.fx]
        self.fx = [i for i in self.fx if check_file(i, emptycheck=True)]
        if len(self.fx) == 0:
            raise ValueError('fx, not found')


    def extract_u1a10(self, fx, outdir=None, gzipped=True, remove=False):
        """Split fx into four groups:
        U1_A10: U1 and A10
        U1_B10: U1 and not A10, (B=C,G,T, not A)
        V1_A10: not U1 and A10, (V=A,C,G, not U)
        V1_B10: not U1 and not A10
        see: IUPAC_code, http://www.bioinformatics.org/sms/iupac.html
        """
        # fname = fq_name(fx)
        x_name = get_fx_name(x, fix_unmap=False)
        if not check_file(fx, emptycheck=True):
            log.error('extract_u1a10() skipped, file not exists, or empty: {}'.format(fx))
            return None
        ftype = Fastx(fx).format # fasta/fastq
        # output dir
        if outdir is None:
            outdir = os.path.join(os.path.dirname(fx), 'U1_A10')
        check_path(outdir)
        # output files
        f_ext = ftype + '.gz' if gzipped else ftype # fastq/fasta
        fout_name = ['U1_A10', 'U1_B10', 'V1_A10', 'V1_B10']
        fout_list = [os.path.join(outdir, '{}.{}.{}'.format(fname, i, f_ext)) \
            for i in fout_name]
        if all(file_exists(fout_list)):
            log.info('FxU1A10() skipped, file exists: {}'.format(fx))
            return fout_list
        # determine output file, writer
        reader = pyfastx.Fastx(fx)
        writer = [xopen(i, 'wt') for i in fout_list]
        def get_writer(u1, a10):
            b1 = '0' if u1 in 'TU' else '1' # U1 
            b2 = '0' if a10 in 'A' else '1' # A10
            return writer[int('0b'+b1+b2, 2)] # binary to int
        if ftype == 'fastq':
            for name,seq,qual,comment in reader:
                if comment:
                    name = name + ' ' + comment
                s = '@{}\n{}\n+\n{}'.format(name, seq, qual)
                w = get_writer(seq[0], seq[9])
                w.write(s+'\n')          
        elif ftype == 'fasta':
            for name,seq,comment in reader:
                if comment:
                    name = name + ' ' + comment
                s = '>{}\n{}'.format(name, seq)
                w = get_writer(seq[0], seq[9])
                w.write(s+'\n')
        else:
            pass
        [w.close() for w in writer]
        # remove 1U 10A files
        if remove:
            file_remove(fout_list, ask=False)
        return fout_list
    
    
    def stat(self):
        """Number of reads
        fx stat
        """
        pass
    
                
    def run_single(self, i):
        """Split fx, by 1U, 10A, single 
        input: index
        """
        return self.extract_u1a10(self.fx[i])
        
        
    def run(self):
        """Check the u1a10 content for fastq files, list
        """
        if len(self.fx) > 1 and self.parallel_jobs > 1:
            with Pool(processes=self.parallel_jobs) as pool:
                pool.map(self.run_single, range(len(self.fx)))
        else:
            for i in range(len(self.fx)):
                self.run_single(i)
        
        # self.run_u1a10(self.fx)

        
        

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
    # input
    if not check_file(fx_in, emptycheck=True):
        raise ValueError('fx_in not exists, {}'.format(fx_in))        
    # force output: fq
    if fx_name.endswith('.fq') or fx_name.endswith('.fastq'):
        fq_out = True
    else:
        fq_out = False
    if file_exists(fx_out):
        log.info('collapse_fx() skipped, file exists: {}'.format(fx_out))
    else:
        d = {}
        fr = pyfastx.Fastx(fx_in)
        for name,seq,qual,comment in fr:
            d[seq] = d.get(seq, 0) + 1
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


def split_fx(fx, outdir=None, min=23, max=29, gzipped=True):
    """Split fastq file by the size 
    23-29 nt
    """
    # f_name = fq_name(fx)
    x_name = get_fx_name(x, fix_unmap=False)
    f_type = Fastx(fx).format
    if outdir is None:
        outdir = os.path.join(os.path.dirname(fx), 'size_select')
    # outfiles
    dir_in = os.path.join(outdir, 'size_select')
    dir_ex = os.path.join(outdir, 'size_exclude')
    fx_name = f_name + '.' + f_type
    if gzipped:
        fx_name += '.gz'
    fx_in = os.path.join(dir_in, fx_name)
    fx_ex = os.path.join(dir_ex, fx_name)
    check_path([dir_in, dir_ex])
    # run
    if not check_file(fx, emptycheck=True):
        log.error('file is empty, split_fx() skipped: {}'.format(fx))
    if all(file_exists([fx_in, fx_ex])):
        log.info('file exists, split_fx() skipped: {}'.format(fx))
    else:
        fr = pyfastx.Fastx(fx)
        with xopen(fx, 'rt') as r, xopen(fx_in, 'wt') as w1, \
            xopen(fx_ex, 'wt') as w2:
            for name,seq,qual,comment in fr:
                if comment:
                    name = name + ' ' + comment
                if f_type == 'fasta':
                    s = '>{}\n{}\n'.format(name, seq)
                elif f_type == 'fastq':
                    s = '@{}\n{}\n+\n{}\n'.format(name, seq, qual)
                else:
                    continue #
                w = w1 if len(seq) in range(min, max+1) else w2
                w.write(s)
    return (fx_in, fx_ex)


def overlap_fx(query, subject, outdir=None, threads=4):
    """Check overlap between fastq files, by sequence
    using command tool: 
    1. seqkit common -s small.fq.gz big.fq.gz
    2. seqkit grep -s -f <(seqkit seq -s small.fq.gz) big.fq.gz # by seq
    """
    if not isinstance(outdir, str):
        outdir = os.path.dirname(query)
    if not check_file([query, subject], emptycheck=True):
        log.error('file not exists, or empty: {}, {}'.format(query, subject))
    else:
        check_path(outdir)
        out_fq = os.path.join(outdir, os.path.basename(query))
        out_log = os.path.join(outdir, 'cmd.log')
        out_cmd = os.path.join(outdir, 'cmd.sh')
        ## slower than seqkit common
        # cmd = ' '.join([
        #     'seqkit grep -s -i -j {}'.format(threads),
        #     '-f <(seqkit seq -s {})'.format(query),
        #     '{}'.format(subject),
        #     '2>{}'.format(out_log),
        #     '| pigz -p {} >{}'.format(threads, out_fq),
        # ])
        cmd = ' '.join([
            'seqkit common -s -j {}'.format(threads),
            '{} {}'.format(query, subject),
            '2> {}'.format(out_log),
            '| pigz -p {} > {}'.format(threads, out_fq)
        ])
        with open(out_cmd, 'wt') as w:
            w.write(cmd + '\n')

        if file_exists(out_fq):
            log.info('overlap_fq() skipped, file exists: {}'.format(out_fq))
        else:
            # os.system(cmd)
            run_shell_cmd(cmd)
        # overlap
        query_c = Fastx(query).number_of_seq()
        ov_c = Fastx(out_fq).number_of_seq()
        msg = 'query: {}, overlap: {}, percent: {:.2f}%'.format(
            query_c, ov_c, ov_c/query_c*100
        )
        print(msg)
        return out_fq
