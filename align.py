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
import pysam
import pandas as pd
import numpy as np
from xopen import xopen
from hiseq.utils.helper import * 
from hiseq.utils.seq import Fastx
from utils import get_fx_name


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('INFO')




class Align(object):
    """Align reads to index, using bowtie
    Align reads to index
    unique + multi (-k 1)
    unique (-m 1)
    multi (-k 2, ids->bam)    
    
    args:
    unique:  -m 1 -v 2 --best -S --un un.fq
    both:    -k 1 -v 2 --best -S --un un.fq
    multi:   -k 2 -v 2 --best -S --un un.fq
    
    Usage:
    Align(fx=in.fq, index=s, outdir='outdir', unique='unique').run()
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()
        
        
    def init_args(self):
        args_init = {
            'fx': None,
            'index': None,
            'outdir': None,
            'unique': 'unique', # unique, multi, both
            'threads': 4,
            'gzipped': True,
            'rm_tmp': False,
        }
        self = update_obj(self, args_init, force=False)
        if not isinstance(self.fx, str):
            raise ValueError('fx expect str, not {}'.format(
                type(self.fx).__name__))
        if not check_file(self.fx, emptycheck=True):
            raise ValueError('fx not exists, or empty: {}'.format(self.fx))
        if not isinstance(self.index, str):
            raise ValueError('index expect str, not {}'.format(
                type(self.index).__name__))
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        check_path(self.outdir)
        self.fx = file_abspath(self.fx)
        self.outdir = file_abspath(self.outdir)
        # output
        self.f_type = Fastx(self.fx).format
        self.f_name = get_fx_name(self.fx, fix_unmap=True)
                
        
    def align(self, fx, outdir=None, fix_outdir=False):
        """Align reads to index
        unique:
        multi:
        both:        
        
        args:
        unique:  -m 1 -v 2 --best -S --un un.fq
        both:    -k 1 -v 2 --best -S --un un.fq
        multi:   -k 2 -v 2 --best -S --un un.fq
        """
        # arguments
        para = '-f' if self.f_type == 'fasta' else '-q'
        para += ' -v 2 --best -S -p {}'.format(self.threads)
        if self.unique == 'unique':
            para += ' -m 1'
        elif self.unique == 'multi':
            para += ' -k 2'
        elif self.unique == 'both':
            para += ' -k 1'
        else:
            pass                
        # output files
        if not isinstance(outdir, str):
            outdir = self.outdir
        if self.unique == 'multi' and not fix_outdir:
            outdir = os.path.join(outdir, 'tmp')
        check_path(outdir)
        prefix = os.path.join(outdir, '{}').format(self.f_name)
        bam = prefix + '.bam'
        bam_log = prefix + '.bowtie.log'
        fx_align = prefix + '.' + self.f_type
        fx_unal = prefix + '.unmap.' + self.f_type
        # cmd
        cmd = ' '.join([
            'bowtie {}'.format(para),
            '--un {}'.format(fx_unal),
            '-x {} {}'.format(self.index, fx),
            '2> {}'.format(bam_log),
            '| samtools view -bhS -F 0x4 -',
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
        if self.gzipped:
            fx_align_gz = fx_align + '.gz'
            fx_unal_gz = fx_unal + '.gz'
            if file_exists(fx_unal):
                gzip_cmd(fx_unal, fx_unal_gz, decompress=False)
            # convert to fq
            self.bam_to_fq(bam, fx_align_gz)
            fout = [fx_align_gz, fx_unal_gz]
        else:
            bam_to_fq(bam, fx_align)
            fout = [fx_align, fx_unal]
        # output
        fout.insert(0, bam)
        return fout
        # return [bam, fx_align_out, fx_unal_out] # return gzipped files
        

    def extract_multi(self, bam, outdir=None):
        """Extract multi aligned reads
        -k 2, 
        ids > 1

        # extract seq from fx, by id
        samtools view -F 0x4 in.bam | cut -f 1 | sort | uniq -d > id_multi
        seqkit grep -n -f id_multi fx > fx_multi
        seqkit grep -v -n -f id_multi fx > fx_unal # unique
        """
        bname = os.path.basename(os.path.splitext(bam)[0])
        if not isinstance(outdir, str):
            outdir = os.path.dirname(bam)
        # output
        prefix = os.path.join(outdir, bname)
        id_multi = prefix + '.multi.id.txt'
        fx_multi = prefix + '.fastq'
        fx_unal = prefix + '.unmap.fastq'
        # get multi ids
        cmd = ' '.join([
            'samtools view -F 0x4 {}'.format(bam),
            '| cut -f 1 | sort | uniq -d > {}'.format(id_multi),
            '&& seqkit grep -n -f {} {} > {}'.format(id_multi, self.fx, fx_multi),
            '&& seqkit grep -v -n -f {} {} > {}'.format(id_multi, self.fx, fx_unal),
        ])
        # save
        cmd_sh = os.path.join(self.outdir, 'cmd_multi.sh')
        with open(cmd_sh, 'wt') as w:
            w.write(cmd + '\n')
        # run
        fx_multi_gz = fx_multi + '.gz'
        fx_unal_gz = fx_unal + '.gz'
        if file_exists(fx_multi) or file_exists(fx_multi_gz):
            log.info('file exists, align skipped: {}'.format(fx_multi))
        else:
            os.system(cmd)
        # gzip output
        if self.gzipped:
            if file_exists(fx_multi):
                gzip_cmd(fx_multi, fx_multi_gz, decompress=False)
            if file_exists(fx_unal):
                gzip_cmd(fx_unal, fx_unal_gz, decompress=False)
            fx_out = [fx_multi_gz, fx_unal_gz]
        else:
            fx_out = [fx_multi, fx_unal]
        return fx_out
        
        
    def revComp(self, s):
        d = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
        s = [d[c] for c in s]
        return ''.join(s[::-1])
        
        
    def bam_to_fq(self, bam, fq=None):
        """
        Convert bam to fastq
        using samtools fastq 
        Note:
        multiple mapping reads
        see also: https://github.com/dpryan79/Answers/blob/master/bioinfoSE_2149/convert.py
        """
        fname = os.path.basename(os.path.splitext(bam)[0])
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
                    s = self.revComp(s)
                    q = q[::-1]
                if fq_format == 'fasta':
                    seq = '>{}\n{}\n'.format(read.query_name, s)
                else:
                    seq = '@{}\n{}\n+\n{}\n'.format(read.query_name, s, q)
                of.write(seq)
            of.close()


    def run(self):
        """Run alignment
        unique: -m 1
        both:   -k 1
        multi:  -k 2 , ids > 1
        """
        fout = self.align(self.fx)
        if self.unique in ['multi']:
            m = self.extract_multi(fout[0], self.outdir)
            # align multi to bam, unique
            self.unique = 'both'
            fout = self.align(fx=m[0], outdir=self.outdir, fix_outdir=True)
        return fout # bam, fq_align, fq_unal


