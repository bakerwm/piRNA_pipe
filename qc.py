#!/usr/bin/env python3



import os
import sys
import logging
import pandas as pd
from multiprocessing import Pool
from hiseq.utils.helper import listfile, update_obj, file_exists, file_copy, Config, run_shell_cmd, file_abspath
from fastx import FxU1A10, FxCount, FxFragSize

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    stream=sys.stdout)
log = logging.getLogger(__name__)
log.setLevel('WARNING')


class PiRNApipeStat(object):
    """Perform the extra analysis
    Parameters:
      x, directory of PiRNApipe() project

    1. count (fx.toml)
    2. U1A10 
    3. length distribution

    overall stat
    """
    def __init__(self, x, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.x = x # project dir
        self.init_args()

    
    def init_args(self):
        args_init = {
            'parallel_jobs': 1,
        }
        self = update_obj(self, args_init, force=False)
        self.x = file_abspath(self.x)
        # config
        config = os.path.join(self.x, 'config', 'config.toml')
        if not file_exists(config):
            raise ValueError('x, missing config: {}'.format(self.x))
        self.prj = Config().load(config) #


    def list_dir(self, u1a10=False, group=False):
        """List all directories
        level-1: 00.raw_data, 01.clean_data, ...
        level-2: /both, /unique, /multi
        include U1A10 or not?!
        
        !!! include 00. -> 12. dirs, 
        !!! exclude config, overlap., ...
        """
        # level-1 00.raw_dat to 09.unmap
        dir1 = listfile(self.x, include_dir=True, recursive=False)
        dir1 = [d for d in dir1 if os.path.basename(d)[0] in '012'] #!!!
        # level-2, both, unique, multi
        if group:
            dir1 += [i + '/both' for i in dir1] + \
                [i + '/unique' for i in dir1] + \
                [i + '/multi' for i in dir1]
        if u1a10:
            dir1ua = [os.path.join(i, 'U1_A10') for i in dir1]
        else:
            dir1ua = [] # empty    
        dirs = dir1 + dir1ua
        return [d for d in dirs if os.path.exists(d)]


    def list_fx(self, u1a10=False, group=False):
        """Return the 00.raw_data to 09.unmap fastq files
        gzipped fq files
        """
        fx_list = []
        dirs = self.list_dir(u1a10, group)
        for d in dirs:
            fx = listfile(d, "*.gz") # recursive ?
            if len(fx) > 0:
                fx_list.extend(fx)
        return fx_list
    
    
    def list_bam(self, u1a10=False, group=True):
        """Return the 04.smRNA to 08.genome bam files
        align, bam files
        
        group: True
        """
        bam_list = []
        dirs = self.list_dir(u1a10, group=True) # level-1 dirs
        for d in dirs:
            bam = listfile(d, "*.bam")
            if len(bam) > 0:
                bam_list.extend(bam)
        return bam_list
    
    
    def list_count_toml(self, d, u1a10=False, group=False):
        """Return the 00.raw_data to 09.unmap *.fx_stat.toml files
        align, bam files
        """
        return listfile(d, "*.fx_stat.toml") # recursive ?
    

    def load_toml(self, x):
        """Convert toml to pd.DataFrame: 
        'sample', 'num_reads', 'fx_type', 'num_seqs'
        
        dict -> pd
        num_reads
        num_seqs
        """
        df = None
        if isinstance(x, str):
            if x.endswith(".toml"):
                d = Config().load(x)
                if len(d) > 0:
                    df = pd.DataFrame.from_dict(d, 'index')
            else:
                log.error('x is not .toml file')
        else:
            log.error('x is not file')
        return df.T # 
    
    
    def stat_dir(self, d):
        """Create summary for level-1 directory
        number of reads, for each category
        total: 
        U1_A10: 
        """
        h = ['sample', 'num_reads', 'fx_type', 'num_seqs']
        dfx = pd.DataFrame(columns=h)
        # level-1: dirs
        group = os.path.basename(d)
        t1 = self.list_count_toml(d)
        if len(t1) == 1:
            df1 = self.load_toml(t1[0])
        else:
            df1 = dfx
        df1.columns = h
        df1['u1a10'] = 'all'
        # level-2: U1A10
        t2 = self.list_count_toml(d+'/U1_A10')
        frames = [self.load_toml(i) for i in t2]
        if len(frames) > 0:
            df2 = pd.concat(frames, axis=0)
            df2.columns = ['sample', 'num_reads', 'fx_type', 'num_seqs']
            df2[['sample', 'u1a10']] = df2['sample'].str.split('.', 1, expand=True)
        else:
            df2 = dfx
            df2['u1a10'] = 'all'
        # combine
        df = pd.concat([df1, df2], 'index')
        df['group'] = group
        # stat
        if len(t1) > 0:
            d_csv = os.path.join(d, 'fx_stat.csv')
            df.to_csv(d_csv, index=False)
        return df
        
    
    def stat_dirs(self):
        """Stat for all directories
        dirs
        """
        dirs = self.list_dir(u1a10=False, group=False)
        frames = [self.stat_dir(i) for i in dirs]
        df = pd.concat(frames, axis=0)
        # save table to 11.stat
        # num_reads
        s1 = os.path.join(self.prj.get('stat_dir'), 'fx_stat.reads.csv')
        df1 = df.drop(['num_seqs', 'fx_type'], axis=1).astype({
            'num_reads': 'int32'})
        df1 = df1.pivot_table(columns='u1a10', 
                              values='num_reads',
                              index=['sample', 'group'])
        df1.reset_index(inplace=True)
        df1 = df1.fillna(0)
        df1.to_csv(s1, index=False)       
        
        # seqs
        s2 = os.path.join(self.prj.get('stat_dir'), 'fx_stat.seqs.csv')
        df2 = df.drop(['num_reads', 'fx_type'], axis=1).astype({
            'num_seqs': 'int32'})
        df2 = df2.pivot_table(columns='u1a10', 
                              values='num_seqs',
                              index=['sample', 'group'])
        df2.reset_index(inplace=True)
        df2 = df2.fillna(0)
        df2.to_csv(s2, index=False)    
        return df
    
    
    def run_stat(self):
        """all analysis
        1.u1a10
        2.count
        3.fragsize
        """
        # 1. U1A10
        bam_list = self.list_bam(u1a10=False) # split bam (group=True)
        FxU1A10(bam_list, parallel_jobs=self.parallel_jobs).run()
        fx_list = self.list_fx(u1a10=False, group=False) # split fx (level-1)
        FxU1A10(fx_list, parallel_jobs=self.parallel_jobs).run()
        
        # 2. count
        fx_list = self.list_fx(u1a10=True, group=False) # count fx
        FxCount(fx_list, parallel_jobs=self.parallel_jobs).run()
        
        # 3. fragsize
        bam_list = self.list_bam(u1a10=True) # size bam
        FxFragSize(bam_list, parallel_jobs=self.parallel_jobs).run()
        fx_list = self.list_fx(u1a10=True, group=False) # size fx
        FxFragSize(fx_list, parallel_jobs=self.parallel_jobs).run()
    
    
    def report(self):
        """Generate report"""
        pkg_dir = os.path.dirname(os.path.realpath(__file__))
        qc_report_r = os.path.join(pkg_dir, 'report.R')
        qc_report_rmd = os.path.join(pkg_dir, 'report.Rmd')        
        # output
        report_dir = self.prj.get('report_dir', None)
        report_rmd = os.path.join(report_dir, os.path.basename(qc_report_rmd))
        report_html = os.path.join(report_dir, 'smRNA_report.html')
        log_out = os.path.join(report_dir, 'log.stdout')
        log_err = os.path.join(report_dir, 'log.stderr')
        # run
        # args: prj_dir, outhtml, template
        cmd = ' '.join([
            'Rscript', qc_report_r, 
            self.x, report_html, report_rmd, 
            '1>', log_out, '2>', log_err])
        cmd_sh = os.path.join(report_dir, 'cmd.sh')
        with open(cmd_sh, 'wt') as w:
            w.write(cmd+'\n')
        if file_exists(report_html):
            log.info('report() skipped, file exists: {}'.format(report_html))
        else:
            file_copy(qc_report_rmd, report_rmd) # copy
            run_shell_cmd(cmd)


    def run(self):
        """Run all stat
        U1A10
        count
        fragsize
        """
        self.run_stat()
        self.stat_dirs()
        self.report()

