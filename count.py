#!/usr/bin/env python3

"""Count reads on features

option-1, chr-level: 
    pysam.idxstats(in.bam)
  
option-2, gene-level:
    fc_anti = hiseq.utils.helper.FeatureCounts(**args_anti)
    fc_anti.run()

## count reads for dirs
input:
{both|unique|multi}/*.bam

output:
{both|unique|multi}/*.count.txt


compatiable with PiRNApipe()
TE, piRC, genome
"""


import os
import sys
import logging
import pybedtools
import pysam
import pandas as pd
import numpy as np
# from piRNA_pipe import PiRNApipe
from hiseq.utils.helper import * # FeatureCounts
from hiseq.utils.bam import Bam


class CountDir(object):
    def __init__(self, prj_dir, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.prj_dir = prj_dir
        self.init_args()
        
    
    def init_args(self):
        args_init = {
            'pipe_obj': None, # PiRNApipe()
            'prj_dir': None, # path to outdir of PiRNApipe()
        }
        self = update_obj(self, args_init, force=False)
#         if isinstance(self.pipe_obj, PiRNApipe):
#             # TE*/piRC*/genome
#             self.te_dir = self.pipe_obj.te_dir
#             self.piRC_dir = self.pipe_obj.piRC_dir
#         elif os.path.exists(self.prj_dir):
        if os.path.exists(self.prj_dir):
            # TE*/piRC*/genome
            self.te_dir = os.path.join(self.prj_dir, '07.te')
            self.piRC_dir = os.path.join(self.prj_dir, '08.piRNA_cluster')
        else:
            raise Exception('missing args: pipe_obj=, or prj_dir=')
        if not file_exists([self.te_dir, self.piRC_dir]):
            raise Exception('07.te_dir, 08.piRNA_cluster, not exists')
            
    
    def count_v1(self, d, group):
        """Count reads for TE/piRC
        chr-level
        pysam.idxstats()
        """
        bam_dir = os.path.join(d, group)
        bam_list = listfile(bam_dir, '*.bam', recursive=False)
        if len(bam_list) == 1:
            bam = bam_list.pop()
            txt = os.path.splitext(bam)[0] + '.count.txt'
            if file_exists(txt):
                log.info('count_v1() skipped, file exists: {}'.format(txt))
            # index
            try:
                Bam(bam).index()
                with open(txt, 'wt') as w:
                    w.write(pysam.idxstats(bam))
            except:
                log.error('count_v1() failed, {}'.format(bam))
        else:
            log.warning('1 bam file supported, {} found'.format(len(bam_list)))


    def count_v2(self, d, group):
        """Count reads for TE/piRC
        gene-level
        FeatureCounts()
        """
        pass
        
        
    def run(self):
        # TE, piRC
        # unique, multi, both
        [self.count_v1(self.te_dir, d) for d in ['unique', 'multi', 'both']]
        [self.count_v1(self.piRC_dir, d) for d in ['unique', 'multi', 'both']]
        
        
        
        
        











