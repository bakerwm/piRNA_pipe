#!/usr/bin/env python3

import os
import sys
import re
import pyfastx
import pandas as pd
import numpy as np
import argparse
from xopen import xopen
from multiprocessing import Pool
from hiseq.fragsize.fragsize import BamFragSizeR1
from hiseq.utils.fastx import Fastx
from hiseq.utils.utils import log, update_obj, run_shell_cmd
from hiseq.utils.file import fx_name, check_path, file_abspath


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
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fh: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fh: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

                
def parse_fx_id(x):
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
        out = [parse_fx_id(i) for i in x]
    else:
        out = None
    return out
                

def cal_freq(x, strandness=True):
    """Calculate the frequency of list
    ['length', 'strand']
    return dataframe
    """
    header = ['length', 'count']
    if isinstance(x, list):
        df = pd.DataFrame(x, columns=header).groupby(['length']).sum().reset_index()
    else:
        df = pd.DataFrame(columns=header)
    return df


def fx_fragsize(x, csv_file=None,  n_max=0):
    """Calculate the length distribution of the fx
    fq
    fa
    >fragsize.csv
    length count
    """
    chunk = 1000000
    frag_size_list = []
    counter = 0
    frames = []    
    with xopen(x) as r:
        for name, seq, qual in readfq(r):
            counter += 1
            n_count = parse_fx_id(name)
            frag_size_list.append([len(seq), n_count])
            # last record
            if n_max > 0 and counter >= n_max:
                log.info('Stopped at limit: {}'.format(n_max))
                break # stop
            # chunk
            if counter > 0 and counter%chunk == 0:
                frames.append(cal_freq(frag_size_list))
                frag_size_list = [] # empty
                log.info('{} : {}'.format('Processed', counter))
        # last chunk
        if len(frag_size_list) > 0:
            frames.append(cal_freq(frag_size_list))
            frag_size_list = [] # empty
            log.info('{} : {}'.format('Processed', counter))
    # overall
    df = pd.concat(frames, axis=0).groupby(['length']).sum().reset_index()
    df['id'] = fx_name(x)
    # save to file
    if isinstance(csv_file, str):
        csv_file = file_abspath(csv_file)
        if os.path.exists(os.path.dirname(csv_file)):
            df.to_csv(csv_file, index=False)
    # output
    return df


def get_args():
    parser = argparse.ArgumentParser(description='hiseq fragsize.py [-o out.csv] <in.fq>')
    parser.add_argument('-o', '--csv-file', dest='csv_file', help='csv file')
    parser.add_argument('fx', help='fastx files')
    return parser


def main():
    args = get_args().parse_args()
    df = fx_fragsize(args.fx, args.csv_file)
    print(df)

if __name__ == '__main__':
    main()

#