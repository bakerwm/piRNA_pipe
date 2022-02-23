#!/usr/bin/env python3
#-*- coding:utf-8 -*-

"""
Extract fastx count

pass read id for count;

input: fastq, fasta
output: fragsize.csv (length,count,id)
"""

import os
import sys
import re
import argparse
from pyfastx import Fastx
from hiseq.utils.utils import update_obj, log, get_date, run_shell_cmd, Config
from hiseq.utils.file import file_exists, file_prefix, file_abspath, check_path


class FxCount(object):
    """Calculate the read count"""
    def __init__(self, fx, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.fx = fx
        if not isinstance(self.fx, str):
            raise ValueError('str expected, {} got'.format(
                type(self.fx).__name__))
        self.fx_name = file_prefix(self.fx)

    
    def is_fa(self, x):
        """Check input is fasta"""
        if x is None:
            x = self.fx
        out = False
        p = re.compile('.f(ast)?a(.gz)?$', flags=re.IGNORECASE)
        # p = re.compile('.f(ast)?q(.gz)?$', flags=re.IGNORECASE)
        if isinstance(x, str):
            m = p.search(x)
            out = isinstance(m, re.Match)
        return out
    
    
    def is_fq(self, x):
        """Check input is fasta"""
        if x is None:
            x = self.fx
        out = False
        # p = re.compile('.f(ast)?a(.gz)?$', flags=re.IGNORECASE)
        p = re.compile('.f(ast)?q(.gz)?$', flags=re.IGNORECASE)
        if isinstance(x, str):
            m = p.search(x)
            out = isinstance(m, re.Match)
        return out

    
    def parse_fx_id(self, x):
        """Extract read count from read id
        if read id in this format:
        1-123 => 123 # count is 123
        abc   => 1 # count 1
        """
        if isinstance(x, str):
            m = re.match('^(\d+)\-(\d+)$', x)
            out = m.group(2) if isinstance(m, re.Match) else 1
            out = int(out) # to int
        elif isinstance(x, list):
            out = [self.parse_fx_id(i) for i in x]
        else:
            out = None
        return out
    
    
    def count_fx(self, x):
        """count reads"""
        if x is None:
            x = self.fx
        u = 0 # unique
        t = 0 # total
        try:
            if self.is_fa(x):
                for name,_,_ in Fastx(x):
                    u += 1
                    t += self.parse_fx_id(name)
            elif self.is_fq(x):
                for name,_,_,_ in Fastx(x):
                    u += 1
                    t += self.parse_fx_id(name)
        except:
            log.error('failed to read fx file: {}'.format(x))
        return [u, t]

    
    def run(self):
        # return '{}\t{}'.format(self.count_fx(self.fx), self.fx_name)
        return [self.fx_name] + self.count_fx(self.fx)

        
def get_args():
    parser = argparse.ArgumentParser(description='hiseq fx_count.py [in.fq (...)]')
    parser.add_argument('fx', nargs='+', help='fastx files')
    return parser


def main():
    args = get_args().parse_args()
    out = []
    for i in args.fx:
        s = FxCount(i).run()
        s = list(map(str, s))
        # print(s)
        out.append('\t'.join(s))
    print('\n'.join(out))


if __name__ == '__main__':
    main()

# EOF