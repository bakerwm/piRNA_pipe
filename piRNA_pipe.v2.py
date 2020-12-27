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
from hiseq.trim.trimmer import Trim
from piRNA_pipe_utils import *


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
    parser.add_argument('-i', '--fq', required=True,
        help='fastq file') 
    parser.add_argument('-o', '--outdir', required=True,
        help='directory to save results')
    parser.add_argument('-w', '--workflow', type=int, default=1,
        help='workflow, 1:4, 1:TE->piRC->Genome; 2:TE->Genome; \
        3:piRC->TE->Genome; 4:piRC->Genome;  default: [1]')
    parser.add_argument('-c', '--collapse', action='store_true',
        help='collapse fastq reads')
    parser.add_argument('-t', '--trimmed', action='store_true',
        help='Input file was clean data, no need to trim adapters')
    parser.add_argument('--subject', default=None,
        help='The target file for overlap')
    parser.add_argument('-p', '--threads', type=int, default=4,
        help='Number of threads to run')
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


class pipeConfig(object):
    """Prepare data for piRNA analysis

    1. index: structureRNA, miRNA, TE, piRNA_cluster, genome, ...
    2. outputdir
    3. pipeline version
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = {
            'outdir': None,
            'smp_name': None,
            'genome': 'dm6',
            'threads': 4,
            'collapse': True,
            'workflow': 1,
            'trimmed': False
        }
        self = update_obj(self, args_local, force=False)
        if self.outdir == None:
            self.outdir = str(pathlib.Path.cwd())
        check_path(self.outdir)

        if self.smp_name == None:
            self.smp_name = fq_name(self.fq, pe_fix=True)

        self.outdir = os.path.join(self.outdir, self.smp_name)

        # update files
        self.init_dirs()
        self.init_files()


    def init_files(self):
        index_dir = ''.join([
            '/home/wangming/data/genome/',
            '{}'.format(self.genome),
            '/bowtie_index'
            ])
        args_f = {
            'index_dir': index_dir,
            'smRNA_index': index_dir + '/smRNA',
            'miRNA_index': index_dir + '/hairpin',
            'te_index': index_dir + '/te',
            'piRC_index': index_dir + '/piRNA_cluster',
            'genome_index': index_dir + '/' + self.genome,
            'config_toml': self.config_dir + '/config.toml'
        }
        self = update_obj(self, args_f, force=True)


    def init_dirs(self):
        """The directory structure
        00.total
        01.collapse
        02.smRNA
        03.miRNA
        04.size_select
        05.TE (could be: piRNA_cluster, ...)
        06.genome
        07.unmap
        08.stat
        """
        arg_dirs = {
            'config_dir': self.outdir + 'config',
            'raw_dir': self.outdir + '/00.raw_data', 
            'clean_dir': self.outdir + '/01.clean_data',
            'collapse_dir': self.outdir + '/02.collapse',
            'smRNA_dir': self.outdir + '/03.smRNA',
            'miRNA_dir': self.outdir + '/04.miRNA',
            'size_dir': self.outdir + '/05.size_select',
            'size_ex_dir': self.outdir + '/05.size_exclude',
            'te_dir': self.outdir + '/06.te', # TE/piR_C/genome
            'piRC_dir': self.outdir + '/07.piRNA_cluster', # piR_C (optional)
            'genome_dir': self.outdir + '/08.genome', # genome (optional)
            'unmap_dir': self.outdir + '/09.unmap',
            'stat_dir': self.outdir + '/10.stat',
            'report_dir': self.outdir + '/11.report'
        }
        self = update_obj(self, arg_dirs, force=True)

        check_path([
            self.raw_dir,
            self.clean_dir,
            self.collapse_dir,
            self.smRNA_dir,
            self.miRNA_dir,
            self.size_dir,
            self.size_ex_dir,
            self.te_dir,
            self.piRC_dir,
            self.genome_dir,
            self.unmap_dir,
            self.stat_dir,
            self.report_dir
            ])


class pipe(object):
    """Run piRNA analysis pipeline

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
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        args_local = pipeConfig(**self.__dict__)
        self = update_obj(self, args_local.__dict__, force=True)


    def prep_raw(self, fx):
        """Copy/symlink raw fastq files
        split files into 1U,10A, add stat
        """
        log.info('00.Copy raw data')
        fq_raw = os.path.join(self.raw_dir, os.path.basename(fx))
        file_symlink(fx, fq_raw)
        # # 1U, 10A
        # fx_list = split_fx_1u10a(fq_raw, remove=False)
        # wrap stat
        fx_stat2(self.raw_dir, recursive=False)
        return fq_raw


    def trim(self, fx, trimmed=False):
        """Cut the adapters from the 3' end
        Trim reads:
        Trim(fq1, outdir, fq2, cut_after_trim='9,-6').run()
        if not:
            copy/links
        """
        fx_clean = os.path.join(self.clean_dir, os.path.basename(fx))
        if trimmed:
            file_symlink(fx, fx_clean)
        else:
            args_local = {
                'fq1': fx,
                'outdir': self.clean_dir,
                'library_type': None,
                'len_min': 15,
                'cut_to_length': 0,
                'recursive': False,
                'parallel_jobs': 1
            }
            trim = Trim(**args_local)
            trim.run()
            file_symlink(trim.clean_fq, fx_clean)
        # check output
        if not check_file(fx_clean, emptycheck=True):
            log.error('Stop at trim(), no reads output: {}'.format(
                fx_clean))
            raise Exception('no reads output: trim()')
        # # 1U, 10A
        # split_fx_1u10a(fx_clean, remove=False)
        # wrap stat
        fx_stat2(self.clean_dir, recursive=False)
        return fx_clean


    def run_collapse(self, fx):
        """Collapse fastq files
        save read count in id line
        """
        log.info('01.collapse fastq')
        fq_collapse = os.path.join(self.collapse_dir, os.path.basename(fx))
        if self.collapse:
            fx_collapse(fx, self.collapse_dir)
        else:
            file_symlink(fx, fq_collapse)
        # # 1U, 10A
        # split_fx_1u10a(fq_collapse, remove=False)
        # wrap stat
        fx_stat2(self.collapse_dir, recursive=False)
        return fq_collapse


    def run_smRNA(self, fx):
        """Map reads to structural RNAs, remove map reads
        unique, multiple
        """
        log.info('02.map to small RNA')
        # alignment, unique + multiple, k=1
        fq_not_smRNA = pipe_align(fx, self.smRNA_index, self.smRNA_dir, 
            self.threads, self.genome, remove_1u10a=False)
        # wrap stat
        fx_stat2(self.smRNA_dir, recursive=False)
        # check output
        if not check_file(fq_not_smRNA, emptycheck=True):
            log.error('Stop at run_smRNA(), no reads output: {}'.format(
                fq_not_smRNA))
            raise Exception('no reads output: run_miRNA()')
        return fq_not_smRNA


    def run_miRNA(self, fx):
        """Map reads to miRNA, remove
        unique, multiple
        """
        log.info('03.map to miRNA')
        fq_not_miRNA = pipe_align(fx, self.miRNA_index, self.miRNA_dir, 
            self.threads, self.genome, remove_1u10a=False)
        # wrap stat
        fx_stat2(self.miRNA_dir, recursive=True)
        # check output
        if not check_file(fq_not_miRNA, emptycheck=True):
            log.error('Stop at run_miRNA(), no reads output: {}'.format(
                fq_not_miRNA))
            raise Exception('no reads output: run_miRNA()')
        return fq_not_miRNA


    def run_size_select(self, fx):
        """Split fastq files by size
        23-29: piRNA candidates
        """
        log.info('04.size select')
        fx_in, fx_ex = splilt_fx_size(fx, outdir=self.size_dir, 
            min=23, max=29, gzipped=True)
        fq_size_in = os.path.join(self.size_dir, os.path.basename(fx_in))
        file_symlink(fx_in, fq_size_in)
        # split_fx_1u10a(fq_size_in, remove = False)
        fx_stat2(self.size_dir, recursive=False)
        # self.run_overlap(fq_size_in) # overlap
        # exclude
        fq_size_ex = os.path.join(self.size_ex_dir, os.path.basename(fx_in))
        file_symlink(fx_ex, fq_size_ex)
        # split_fx_1u10a(fq_size_ex, remove=False)
        fx_stat2(self.size_ex_dir, recursive=False)
        # check output
        if not check_file(fq_size_in, emptycheck=True):
            log.error('Stop at run_size_select(), no reads output: {}'.format(
                fq_size_in))
            raise Exception('no reads output: run_size_select()')
        return fq_size_in


    def run_TE(self, fx):
        """Map reads to TE consensus
        """
        log.info('05.map to TE')
        fq_not_te = pipe_align(fx, self.te_index, self.te_dir, 
            self.threads, self.genome, unique_multi='all')
        fq_te = os.path.join(self.te_dir, self.smp_name + '.fastq.gz')
        # self.run_overlap(fq_te)
        # wrap stat
        fx_stat2(self.te_dir, recursive=False)
        # check output
        if not check_file(fq_not_te, emptycheck=True):
            log.error('Stop at run_TE(), no reads output: {}'.format(
                fq_not_te))
            raise Exception('no reads output: run_TE()')
        return fq_not_te


    def run_piRC(self, fx):
        """Map reads to piRNA cluster
        mapping
        """
        log.info('06.map to piRNA cluster')
        fq_not_piRC = pipe_align(fx, self.piRC_index, self.piRC_dir, 
            self.threads, self.genome, unique_multi='all')
        fq_piRC = os.path.join(self.piRC_dir, self.smp_name + '.fastq.gz')
        # self.run_overlap(fq_piRC)
        # wrap stat
        fx_stat2(self.piRC_dir, recursive=False)
        # check output
        if not check_file(fq_not_piRC, emptycheck=True):
            log.error('Stop at run_piRC(), no reads output: {}'.format(
                fq_not_piRC))
            raise Exception('no reads output: run_piRC()')
        return fq_not_piRC


    def run_genome(self, fx):
        """map reads to genome
        non-TE reads to genome
        """
        log.info('07.map to genome')
        fq_not_genome = pipe_align(fx, self.genome_index, self.genome_dir, 
            self.threads, self.genome, unique_multi='all')
        fq_genome = os.path.join(self.genome_dir, self.smp_name + '.fastq.gz')
        # self.run_overlap(fq_genome)
        # wrap stat
        fx_stat2(self.genome_dir, recursive=False)
        # check output
        if not check_file(fq_not_genome, emptycheck=True):
            log.error('Stop at run_genome(), no reads output: {}'.format(
                fq_not_genome))
            raise Exception('no reads output: run_genome()')
        # save unmap
        fx_name = os.path.basename(fx).replace('.unmap', '')
        fx_unmap = os.path.join(self.unmap_dir, fx_name)
        file_symlink(fq_not_genome, fx_unmap)
        # split_fx_1u10a(fx_unmap)
        fx_stat2(self.unmap_dir, recursive=False)
        return fq_not_genome


    ############################
    ## collapse, RNA species
    ############################
    def run_unique_reads(self):
        """Collapse reads in level-1
        count reads
        """
        log.info('collapse reads to piRNA species')
        fx1 = self.get_subfiles(level=1)
        fx2 = self.get_subfiles(level=1, group='overlap')
        fx3 = self.get_subfiles(level=2)
        fx4 = self.get_subfiles(level=2, group='overlap')
        fx_subfiles = fx1 + fx2 + fx3 + fx4
        for fx in fx_subfiles:
            fx_collapse_dir = os.path.join(os.path.dirname(fx), 'collapse')
            fx_collapse(fx, fx_collapse_dir)

        # fx stat
        dir1 = self.get_subdirs(level=1)
        dir2 = self.get_subdirs(level=1, group='overlap')
        dir3 = self.get_subdirs(level=2)
        dir4 = self.get_subdirs(level=2, group='overlap')
        fx_subdirs = dir1 + dir2 + dir3 + dir4
        for fx in fx_subdirs:
            fx_collapse_dir = os.path.join(fx, 'collapse')
            fx_stat2(fx_collapse_dir)


    ############################
    ## statistics
    ## 1. 1U, 10A
    ## 2. overlap with subject 
    ############################
    def get_subdirs(self, level=1, group=None):
        """Return the working directories
        group: overlap, collapse
        level-1: 00.raw_data/
        level-2: 00.raw_data/1U_10A/
        """
        # level-1 00.raw_dat to 09.unmap
        fx_dirs = listfile(self.outdir, include_dir=True)
        # level-2 1u10a dirs
        fx_1u10a = [os.path.join(i, '1U_10A') for i in fx_dirs]
        fx_1u10a = [i for i in fx_1u10a if os.path.exists(i)]
        # output
        out = fx_1u10a if level == 2 else fx_dirs
        # add group
        if isinstance(group, str):
            out = [os.path.join(i, group) for i in out]
        return out


    def get_subfiles(self, level=1, group=None):
        """Return the 00.raw_data to 09.unmap fastq files
        group: overlap, collapse
        level-1: 00.raw_data/*fastq.gz
        level-2: 00.raw_data/1U_10A/*.fastq.gz
        """
        fx_list = []
        fx_dirs = self.get_subdirs(level, group)
        for fx in fx_dirs:
            fx_gz = listfile(fx, '*.gz', recursive=False)
            if len(fx_gz) > 0:
                fx_list.extend(fx_gz)

        return fx_list


    def run_1u10a(self):
        """Check the 1U10A content for the subfiles
        level-1: for 00.raw_data/*.fastq.gz 
        """
        fx_subfiles = self.get_subfiles(level=1)
        for fx in fx_subfiles:
            fx_1u10a = os.path.join(os.path.dirname(fx), '1U_10A')
            split_fx_1u10a(fx, remove=False)
            fx_stat2(fx_1u10a)


    def run_overlap(self):
        """Check fastq files overlap with subject
        level-1: 00.raw_data/*fastq.gz
        level-2: 00.raw_data/1U_10A/*.fastq.gz
        """
        # level-1: single files
        fx_subfiles = self.get_subfiles(level=1)
        for fx in fx_subfiles:
            fx_ov = os.path.join(os.path.dirname(fx), 'overlap')
            fq_overlap(fx, self.subject, fx_ov)
            fx_stat2(fx_ov)

        # level-2: multi files
        fx_1u10a_dirs = self.get_subdirs(level=2)
        fx_1u10a_files = self.get_subfiles(level=2)
        for fx in fx_1u10a_files:            
            fx_ov = os.path.join(os.path.dirname(fx), 'overlap')
            fq_overlap(fx, self.subject, fx_ov)
        for fx in fx_1u10a_dirs:
            fx_ov = os.path.join(fx, 'overlap')
            fx_stat2(fx_ov)


    def read_fx_stat(self, x):
        """Convert toml to pd.DataFrame
        dict -> pd

        num_seqs
        """
        df = None
        if isinstance(x, str):
            if x.endswith(".toml"):
                df = pd.DataFrame.from_dict(Toml(x).to_dict())
            else:
                log.error('x is not .toml file')
        else:
            log.error('x is not file')
        return df


    def read_fx_level1(self):
        """Read the fx stat in 00.raw_data/ level
        count
        group
        """
        fx_dirs = self.get_subdirs(level=1)
        fx_list = [i + '/fx_stat.toml' for i in fx_dirs]
        fx_list = [i for i in fx_list if file_exists(i)]
        fx_frames = [self.read_fx_stat(i) for i in fx_list]
        df = pd.concat(fx_frames, axis=0).reset_index()
        df.columns = ['sample', 'num_seqs']
        g_list = [os.path.basename(os.path.dirname(i)) for i in fx_list]
        df['group'] = g_list
        df['1u10a'] = 'all'
        return df


    def read_fx_level2(self):
        """Read 1U10A in 00.raw_data/1U10A/
        count
        group
        1u10a
        """
        fx_dirs = self.get_subdirs(level=2)
        fx_list = [i + '/fx_stat.toml' for i in fx_dirs]
        fx_list = [i for i in fx_list if file_exists(i)]
        fx_frames = []
        for i in fx_list:
            df = self.read_fx_stat(i)
            df.reset_index(inplace=True)
            df.columns = ['1u10a', 'num_seqs']
            df['group'] = pathlib.Path(i).parts[-3]
            df['sample'] = self.smp_name
            fx_frames.append(df)
        df = pd.concat(fx_frames, axis=0)
        return df


    def read_fx_level1x(self):
        """Read the fx stat in 00.raw_data/overlap level
        count
        group
        """
        fx_dirs = self.get_subdirs(level=1)
        fx_list = [i + '/overlap/fx_stat.toml' for i in fx_dirs]
        fx_list = [i for i in fx_list if file_exists(i)]
        fx_frames = []
        for i in fx_list:
            df = self.read_fx_stat(i)
            df.reset_index(inplace=True)
            df.columns = ['overlap', 'num_seqs']
            df['group'] = pathlib.Path(i).parts[-3]
            df['sample'] = self.smp_name
            df['1u10a'] = 'all'
            fx_frames.append(df)
        df = pd.concat(fx_frames, axis=0)
        return df


    def read_fx_level2x(self):
        """Read 1U10A in 00.raw_data/1U_10A/overlap/
        count
        group
        1u10a
        """
        fx_dirs = self.get_subdirs(level=2)
        fx_list = [i + '/overlap/fx_stat.toml' for i in fx_dirs]
        fx_list = [i for i in fx_list if file_exists(i)]
        fx_frames = []
        for i in fx_list:
            df = self.read_fx_stat(i)
            df.reset_index(inplace=True)
            df.columns = ['overlap', 'num_seqs']
            df['group'] = pathlib.Path(i).parts[-4]
            df['sample'] = self.smp_name
            fx_frames.append(df)
        df = pd.concat(fx_frames, axis=0)
        return df


    def read_fx_collapse(self):
        """Read 00.raw_data/collapse, 
        00.raw_data/1U_10A/collapse
        for reads
        """
        # level-1
        fx_dirs = self.get_subdirs(level=1)
        fx_list = [i + '/collapse/fx_stat.toml' for i in fx_dirs]
        fx_list = [i for i in fx_list if file_exists(i)]
        fx_frames = []
        for i in fx_list:
            df = self.read_fx_stat(i)
            df.reset_index(inplace=True)
            df.columns = ['sample', 'num_seqs']
            df['group'] = pathlib.Path(i).parts[-3]
            df['1u10a'] = 'all'
            fx_frames.append(df)
        df1 = pd.concat(fx_frames, axis=0)

        # level-2
        fx_dirs = self.get_subdirs(level=2)
        fx_list = [i + '/collapse/fx_stat.toml' for i in fx_dirs]
        fx_list = [i for i in fx_list if file_exists(i)]
        fx_frames = []
        for i in fx_list:
            df = self.read_fx_stat(i)
            df.reset_index(inplace=True)
            df.columns = ['unique', 'num_seqs']
            df['group'] = pathlib.Path(i).parts[-4]
            fx_frames.append(df)
        df2 = pd.concat(fx_frames, axis=0)
        df2[['sample', '1u10a']] = df2.unique.str.split('.', 1, expand=True)
        df2 = df2.drop(['unique'], axis=1)

        # combine level1, level2
        df3 = pd.concat([df1, df2], axis=0)
        df3x = df3.pivot_table(columns='1u10a', values='num_seqs', index=['sample', 'group'])
        df3x.reset_index(inplace=True)

        # save to file
        stat_txt = os.path.join(self.stat_dir, 'fx_stat.collapse.txt')
        df3x.to_csv(stat_txt, index=False)


    def read_fx_collapse2(self):
        """Read 00.raw_data/overlap/collapse, 
        00.raw_data/1U_10A/overlap/collapse
        for reads
        """
        # level-1
        fx_dirs = self.get_subdirs(level=1, group='overlap')
        fx_list = [i + '/collapse/fx_stat.toml' for i in fx_dirs]
        fx_list = [i for i in fx_list if file_exists(i)]
        fx_frames = []
        for i in fx_list:
            df = self.read_fx_stat(i)
            df.reset_index(inplace=True)
            df.columns = ['overlap', 'num_seqs']
            df['group'] = pathlib.Path(i).parts[-4]
            df['1u10a'] = 'all'
            fx_frames.append(df)
        df1 = pd.concat(fx_frames, axis=0)        
        df1[['overlap', 'sample']] = df1.overlap.str.split('.', 1, expand=True)

        # level-2
        fx_dirs = self.get_subdirs(level=2, group='overlap')
        fx_list = [i + '/collapse/fx_stat.toml' for i in fx_dirs]
        fx_list = [i for i in fx_list if file_exists(i)]
        fx_frames = []
        for i in fx_list:
            df = self.read_fx_stat(i)
            df.reset_index(inplace=True)
            df.columns = ['overlap', 'num_seqs']
            df['group'] = pathlib.Path(i).parts[-5]
            fx_frames.append(df)
        df2 = pd.concat(fx_frames, axis=0)
        df2[['overlap', 'sample', '1u10a']] = df2.overlap.str.split('.', 2, expand=True)

        # combine level1, level2
        df3 = pd.concat([df1, df2], axis=0)
        df3x = df3.pivot_table(columns='1u10a', values='num_seqs', index=['sample', 'group', 'overlap'])
        df3x.reset_index(inplace=True)

        # save to file
        stat_txt = os.path.join(self.stat_dir, 'fx_stat.overlap.collapse.txt')
        df3x.to_csv(stat_txt, index=False)


    def run_stat(self):
        """Create stat report for the alignment
        1. number of reads in categories
        2. number of reads in overlap
        """
        # 1. reads in each categories
        df1 = self.read_fx_level1()
        df2 = self.read_fx_level2()
        df = pd.concat([df1, df2], axis=0)
        df['1u10a'] = [re.sub(self.smp_name + '.', '', i) for i in  df['1u10a']]
        dfx = df.pivot_table(columns='1u10a', values='num_seqs', index=['sample', 'group'])
        dfx.reset_index(inplace=True)
        f1_stat_txt = os.path.join(self.stat_dir, 'fx_stat.txt')
        dfx.to_csv(f1_stat_txt, index=False)

        # 2. including overlap reads
        df1x = self.read_fx_level1x()
        df2x = self.read_fx_level2x()
        df1x[['overlap', 'sample']] = df1x.overlap.str.split('.', 1, expand=True)
        df2x[['overlap', 'sample', '1u10a']] = df2x.overlap.str.split('.', 2, expand=True)
        df3 = pd.concat([df1x, df2x], axis=0)
        df3x = df3.pivot_table(columns='1u10a', values='num_seqs', index=['sample', 'group', 'overlap'])
        df3x.reset_index(inplace=True)
        f2_stat_txt = os.path.join(self.stat_dir, 'fx_stat.overlap.txt')
        df3x.to_csv(f2_stat_txt, index=False)

        # 3. collapse piRNAs for each categories
        self.read_fx_collapse()

        # 4. collapse piRNAs in overlap
        self.read_fx_collapse2()


    def run(self):
        fq_raw = self.prep_raw(self.fq)
        fq_clean = self.trim(fq_raw, trimmed=True)
        fq_collapse = self.run_collapse(fq_raw)
        fq_not_smRNA = self.run_smRNA(fq_collapse)
        fq_not_miRNA = self.run_miRNA(fq_not_smRNA)
        fq_size = self.run_size_select(fq_not_miRNA)
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

        ## for statistics
        self.run_1u10a()
        self.run_unique_reads()
        self.run_overlap()
        self.run_stat()


def fq_overlap(query, subject, outdir):
    """Compare, overlap between fastq files
    using list, set
    """
    check_path(outdir)
    d = {}
    with xopen(subject, 'rt') as r:
        for n, s, q, m in readfq(r):
            # version-1
            # check: 1-23 nt, for overlap
            s = s[:23]
            d[s] = d.get(s, 0) + 1
    # all unique subject sequences, 1-23 nt
    sub = set(d.keys())

    # output files
    fx_name = fq_name(query)
    fx_in = os.path.join(outdir, 'overlap.'+os.path.basename(query))
    fx_ex = os.path.join(outdir, 'not_overlap.'+os.path.basename(query))
    if all(file_exists([fx_in, fx_ex])):
        log.info('fq_overlap() skipped, file exists: {}'.format(fx_name))
    else:
        with xopen(query, 'rt') as r, \
            xopen(fx_in, 'wt') as w1, \
            xopen(fx_ex, 'wt') as w2:
            for n, s, q, m in readfq(r):
                fq = '\n'.join(['@'+n, s, '+', q])
                s1 = s[:23]
                w = w1 if s1 in sub else w2
                w.write(fq + '\n') 

    return fx_in


def main():
    # args = get_args()
    args = vars(get_args())
    pipe(**args).run()


if __name__ == '__main__':
    main()


## EOF