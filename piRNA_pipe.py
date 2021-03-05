#!/usr/bin/env python3

"""
piRNA analysis (small RNAseq) 

Example:

1. Calculate overlap between files
piRNA_pipe.py -i in.fq -o outdir -g dm6 -t -p 4 -j 4 -s a.fq b.fq c.fq



date: 2020-12-27

1. support collapsed reads, counting, RNA seqs, reads


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
import Levenshtein
from multiprocessing import Pool
from xopen import xopen
from hiseq.trim.trimmer import Trim
from hiseq.fragsize.fragsize import BamFragSize
from hiseq.utils.helper import listfile
from hiseq.utils.seq import Fastx
from piRNA_pipe_utils import *
from align import Align
from fastx import collapse_fx, split_fx, overlap_fx
from utils import get_args, PipeConfig, get_fx_name
from qc import PiRNApipeStat


logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('WARNING')


        
class PiRNApipe(object):
    """Run piRNA analysis pipeline
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
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        args_local = PipeConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)
        Config().to_toml(self.__dict__, self.config_toml)


    def prep_raw(self, fq):
        """Copy/symlink raw fastq files
        split files into 1U,10A, add stat
        """
        log.info('00.Copy raw data')
        file_symlink(fq, self.fq_raw)
        return self.fq_raw


    def run_trim(self, fq, trimmed=True):
        """Cut the adapters from the 3' end
        Trim reads:
        Trim(fq1, prj_dir, fq2, cut_after_trim='9,-6').run()
        if not:
            copy/links
        """
        log.info('01.Trim adapters')
        try:
            if trimmed:
                file_symlink(fq, self.fq_clean)
            else:
                args_local = {
                    'fq1': fq,
                    'outdir': self.clean_dir,
                    'library_type': None,
                    'len_min': 18,
                }
                trim = Trim(**args_local)
                trim.run()
                file_symlink(trim.clean_fq, self.fq_clean)
        except:
            log.error('run_trim() failed')
            fq_clean = None
        return self.fq_clean


    def run_collapse(self, fq):
        """Collapse fastq files! required
        save read count in id line
        """
        log.info('02.Collapse fastq')
        if self.collapse:
            try:
                out_fq = collapse_fx(fq, self.fq_collapse)
            except:
                log.error('run_collapse() failed')
                fq_collapse = None
        else:
            file_symlink(fq, self.fq_collapse)
        return self.fq_collapse


    def run_overlap(self, fq):
        """Check overlap between fq and subject files: could be N files
        Input fastq
        """
        log.info('03.Overlap with other RNAs')
        if self.subject:
            try:
                overlap_fx(fq, self.subject, self.overlap_dir)
            except:
                log.error('run_overlap() failed')
                fq_overlap = None
        else:
            file_symlink(fq, self.fq_overlap)
        return self.fq_overlap


    def run_smRNA(self, fq):
        """Map reads to structural RNAs, remove map reads
        unique + multiple, both
        """
        log.info('04.Map to small RNA')
        fout = self.run_align(fq, 'smRNA', 'both')
        return fout[2] # bam, fq_align, fq_unal


    def run_miRNA(self, fq):
        """Map reads to miRNA, remove
        unique, multiple
        """
        log.info('05.Map to miRNA')
        fout = self.run_align(fq, 'miRNA', 'both')
        return fout[2] # bam, fq_align, fq_unal


    def run_size_select(self, fq):
        """Split fastq files by size
        23-29: piRNA candidates
        """
        log.info('06.Size select')
        try:
            fout = split_fx(fq, self.size_dir, min=23, max=29) # in, ex
            file_symlink(fout[0], self.fq_size_in)
            file_symlink(fout[1], self.fq_size_ex)
        except:
            log.error('run_size_select() failed')
            self.fq_size_in = None
        return self.fq_size_in


    def run_TE(self, fq):
        """Map reads to TE consensus
        unique reads
        """
        log.info('07.map to TE')
        fout = self.run_align(fq, 'te', 'both')
        self.run_align(fq, 'te', 'unique')
        self.run_align(fq, 'te', 'multi')
        return fout[2] # bam, fq_align, fq_unal


    def run_piRC(self, fq):
        """Map reads to piRNA cluster
        mapping
        """
        log.info('08.Map to piRNA cluster')
        fout = self.run_align(fq, 'piRC', 'both')
        self.run_align(fq, 'piRC', 'unique')
        self.run_align(fq, 'piRC', 'multi')
        return fout[2] # bam, fq_align, fq_unal


    def run_genome(self, fq):
        """map reads to genome
        non-TE reads to genome
        """
        log.info('09.Map to genome')
        fout = self.run_align(fq, 'genome', 'both')
        self.run_align(fq, 'genome', 'unique')
        self.run_align(fq, 'genome', 'multi')
        # link unal
        file_symlink(fout[2], self.fq_unal)
        return fout[2] # bam, fq_align, fq_unal


    def run_align(self, fq, group, unique='both'):
        """Align reads to index"""
        args_d = {
            'smRNA': {
                'index': self.smRNA_index, 
                'outdir': os.path.join(self.smRNA_dir, unique),
            },
            'miRNA': {
                'index': self.miRNA_index, 
                'outdir': os.path.join(self.miRNA_dir, unique),
            },
            'te': {
                'index': self.te_index, 
                'outdir': os.path.join(self.te_dir, unique),
            },
            'piRC': {
                'index': self.piRC_index, 
                'outdir': os.path.join(self.piRC_dir, unique),
            },
            'genome': {
                'index': self.genome_index, 
                'outdir': os.path.join(self.genome_dir, unique),
            },
        }
        # check args
        args = args_d.get(group, {})
        args.update({
            'fx': fq,
            'unique': unique, 
            'threads': self.threads,
        })
        try:
            fout = Align(**args).run() # bam, fq_align, fq_unal
            # create symlink, move-up level-1
            fq_dir = os.path.dirname(os.path.dirname(fout[1]))
            fq_align = os.path.join(fq_dir, os.path.basename(fout[1]))
            file_symlink(fout[1], fq_align)
            fout[1] = fq_align # update
        except:
            log.error('run_{}(unique={}) failed'.format(group, unique))
            fout = (None, None, None)
        return fout
        
    
    def run(self):
        # main
        # raw->clean->collapse->overlap->smRNA->miRNA->size_select
        fq_size = self.run_size_select(
            self.run_miRNA(
                self.run_smRNA(
                    self.run_overlap(
                        self.run_collapse(
                            self.run_trim(
                                self.prep_raw(self.fq)))))))
        # switch workflow
        if self.workflow == 1:
            # workflow-1: smRNA->miRNA->size->TE->piRC->genome
            fq_not_genome = self.run_genome(
                self.run_piRC(
                    self.run_TE(fq_size)))
        elif self.workflow == 2:
            # workflow-2: smRNA->miRNA->size->TE->genome
            fq_not_genome = self.run_genome(
                self.run_TE(fq_size))
        elif self.workflow == 3:
            # workflow-1: smRNA->miRNA->size->piRC->TE->genome
            fq_not_genome = self.run_genome(
                self.run_TE(
                    self.run_piRC(fq_size)))
        elif self.workflow == 4:
            # workflow-1: smRNA->miRNA->size->piRC->TE->genome
            fq_not_genome = self.run_genome(
                    self.run_piRC(fq_size))
        else:
            raise Exception('workflow expect [1:4], default: 1')
        # stat
        PiRNApipeStat(self.prj_dir).run()


def overlap_subject(p):
    """Run overlap, with subjects
    p is the object of PiRNApipe
    """
    if not isinstance(p, PiRNApipe):
        log.error('run_overlap() failed, p expect PiRNApipe, got {}'.format(
            type(p).__name__))
    else:
        args_local = p.__dict__
        if len(p.subject_list) > 0:
            for sub in p.subject_list:
                # s_name = fq_name(sub)
                s_name = get_fx_name(sub)
                args_ov = {
                    'fq': p.fq_collapse,
                    'trimmed': True,
                    'collapse': True,
                    'force_overlap': True,
                    'subject': sub,
                    'outdir': os.path.join(p.prj_dir, 'overlap', s_name)
                }
                args_local.update(args_ov)
                pv = PiRNApipe(**args_local)
                pv.run()

        
def main():
    args = get_args()
    # Not for overlap
    args_local = vars(args).copy()
    args_local['force_overlap'] = False
    p = PiRNApipe(**args_local)
    p.run()
    
    # For overlap
    overlap_subject(p)


if __name__ == '__main__':
    main()


## EOF