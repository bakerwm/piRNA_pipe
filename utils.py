#!/usr/bin/env python3
"""Helper functions for piRNA_pipe

get_args(), for arguments 
PipeConfig(), for preparing data, files, dirs, 
"""

import os
import argparse
from hiseq.utils.helper import * # update_obj, fq_name, file_abspath, check_path, check_file






def get_fx_name(x, fix_unmap=True):
    """Parse the name of fastx
    a.fq
    a.fq.gz
    a.unmap.fq
    a.unmap.fq.gz
    """
    if isinstance(x, str):
        x_name = os.path.basename(x)
        # remove .gz
        if x_name[-3:] == '.gz':
            x_name = x_name[:-3]
        # remove ext:
        x_name = os.path.splitext(x_name)[0]
        # remove unmap
        if x_name[-6:] == '.unmap' and fix_unmap:
            x_name = x_name[:-6]
        return x_name
        
        
def rev_comp(s):
    d = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    s = [d[c] for c in s]
    return ''.join(s[::-1])




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
    parser.add_argument('-g', '--genome', default='dm6',
                        help='reference genome, default: [dm6]')
    parser.add_argument('-t', '--trimmed', action='store_true',
                        help='No need to trim adapters')
    parser.add_argument('-s', '--subject-list', nargs='+', dest='subject_list',
                        default=[],
                        help='list of fastq files, for overlap, default: []')
    parser.add_argument('--ov-type', dest='ov_type', type=int, default=1,
                        help='overllap type, 1=perfect match, \
                        2=allow mismatch, default: [1]')
    parser.add_argument('-w', '--workflow', type=int, default=1,
                        help='workflow, 1:4, 1:TE->piRC->Genome; \
                        2:TE->Genome; 3:piRC->TE->Genome; \
                        4:piRC->Genome;  default: [1]')
    parser.add_argument('-p', '--threads', type=int, default=4,
                        help='Number of threads to run')
    parser.add_argument('-j', '--parallel-jobs', type=int, default=1,
                        dest='parallel_jobs',
                        help='number of jobs to run in parallel')
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


class PipeConfig(object):
    """
    Check arguments, default arguments
    1. index: structureRNA, miRNA, TE, piRNA_cluster, genome, ...
    2. outputdir
    3. pipeline version
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_local = {
            'fq': None,
            'outdir': None,
            'smp_name': None,
            'genome': 'dm6',
            'threads': 4,
            'parallel_jobs': 1,
            'collapse': True,
            'workflow': 1,
            'trimmed': True,
            'force_overlap': False,
            'subject': None,
            'subject_list': None,
            'ov_type': 1,
            'genome_path': '~/data/genome'
        }
        self = update_obj(self, args_local, force=False)
        # smp_name, outdir, ...
        if isinstance(self.fq, str):
            if not check_file(self.fq, emptycheck=True):
                raise ValueError('fq not exists, or empty: {}'.format(self.fq))
        else:
            raise ValueError('fq expect str, not {}'.format(
                type(self.fq).__name))
        # smp_name
        if not isinstance(self.smp_name, str):
            self.smp_name = fq_name(self.fq, pe_fix=False)
        self.fx_name = self.smp_name + '.fastq.gz' # force fastq output
        # outdir
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        self.outdir = file_abspath(self.outdir) # absolute path
        self.prj_dir = os.path.join(self.outdir, self.smp_name) # update: outdir !!!        
        # genome_path
        if not isinstance(self.genome_path, str):
            raise ValueError('genome_path expect str, got {}'.format(
                type(self.genome_path).__name__))
        self.genome_path = os.path.expanduser(self.genome_path)
        # extra arguments
        if not isinstance(self.genome, str):
            raise ValueError('genome, expect str, not {}'.format(
                type(self.genome).__name__))
        # threads, parallel_jobs, collapse, trimmed, ...
        self.fq = file_abspath(self.fq)
        self.outdir = file_abspath(self.outdir)
        # update files
        self.init_overlap()
        self.init_dirs()
        self.init_files()
        # log files
        msg = '\n'.join([
            '{:>15s}: {}'.format('fq', self.fq),
            '{:>15s}: {}'.format('outdir', self.outdir),
            '{:>15s}: {}'.format('smp_name', self.smp_name),
            '{:>15s}: {}'.format('genome', self.genome),
            '{:>15s}: {}'.format('subject', self.subject),
            '{:>15s}: {}'.format('overlap', 'yes' if self.force_overlap else 'no'),
        ])
        print('-'*80)
        print(msg)
        # index
        self.init_index()

    
    def init_overlap(self):
        """Check the extra arguments, for overlap
        subject
        overlap
        """        
        # subject, fastq for bowtie index
        if self.subject is None:
            pass
        elif isinstance(self.subject, str):
            if not file_exists(self.subject):
                self.subject = None
        else:
            raise ValueError('subject, expect str, got {}'.format(
                type(self.subject).__name__))
        # subject_list, from args
        if self.subject_list is None:
            self.subject_list = []
        if isinstance(self.subject_list, str):
            self.subject_list = [self.subject_list]
        if not isinstance(self.subject_list, list):
            self.subject_list = [] # force
        self.subject_list = [i for i in self.subject_list if file_exists(i)]            
        # overlap dir
        if self.force_overlap and isinstance(self.subject, str):
            pass
            # s_name = fq_name(self.subject)
            # self.prj_dir = os.path.join(self.prj_dir, 'overlap', s_name) #update
        else:
            self.subject = None # force not subject
        

    def init_index(self):
        """Checkt the bowtie indexes
        genome_path: ~/data/genome (default)
        require: smRNA, miRNA, te, piRC, genome
        """
        index_dir = os.path.join(self.genome_path, self.genome, 'bowtie_index')
        if not os.path.exists(index_dir): 
            raise ValueError('genome_path not exists: {}'.format(index_dir))
        self.index_dir = index_dir
        arg_index = {
            'smRNA_index': index_dir + '/smRNA',
            'miRNA_index': index_dir + '/hairpin',
            'te_index': index_dir + '/te',
            'piRC_index': index_dir + '/piRNA_cluster',
            'genome_index': index_dir + '/' + self.genome
        }
        # check index exists
        fc = []
        print('-'*80)
        print('Check the bowtie indexes:')
        for k, v in arg_index.items():
            f = v + '.1.ebwt'
            ff = 'ok' if os.path.exists(f) else 'failed'
            fc.append(ff == 'ok') # check
            msg = '{:>15s}: {:<6s} {:s}'.format(k, ff, v)
            print(msg)
        print('-'*80)
        if not all(fc):
            raise ValueError('Missing indexes, check above message')
        self = update_obj(self, arg_index, force=True)
        

    def init_dirs(self):
        """The directory structure
        updated: 2021-01-10
        """
        args_dirs = {
            'config_dir': self.prj_dir + '/config',
            'raw_dir': self.prj_dir + '/00.raw_data', 
            'clean_dir': self.prj_dir + '/01.clean_data',
            'collapse_dir': self.prj_dir + '/02.collapse',
            'overlap_dir': self.prj_dir + '/03.overlap',
            'smRNA_dir': self.prj_dir + '/04.smRNA',
            'miRNA_dir': self.prj_dir + '/05.miRNA',
            'size_dir': self.prj_dir + '/06.size_select',
            'size_ex_dir': self.prj_dir + '/06.size_exclude',
            'te_dir': self.prj_dir + '/07.te', # TE/piR_C/genome
            'piRC_dir': self.prj_dir + '/08.piRNA_cluster', # piR_C (optional)
            'genome_dir': self.prj_dir + '/09.genome', # genome (optional)
            'unmap_dir': self.prj_dir + '/10.unmap',
            'stat_dir': self.prj_dir + '/11.stat',
            'report_dir': self.prj_dir + '/12.report'
        }
        self = update_obj(self, args_dirs, force=True)
        check_path(list(args_dirs.values()))
                
            
    def init_files(self):
        """Default files
        config.toml
        ...
        """
        args_files = {
            'config_toml': self.config_dir + '/config.toml',
            'fq_raw': os.path.join(self.raw_dir, self.fx_name),
            'fq_clean': os.path.join(self.clean_dir, self.fx_name),
            'fq_collapse': os.path.join(self.collapse_dir, self.fx_name),
            'fq_overlap': os.path.join(self.overlap_dir, self.fx_name),
            'fq_size_in': os.path.join(self.size_dir, self.fx_name),
            'fq_size_ex': os.path.join(self.size_ex_dir, self.fx_name),
            'fq_unal': os.path.join(self.unmap_dir, self.fx_name),
        }
        self = update_obj(self, args_files, force=True)
