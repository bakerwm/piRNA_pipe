#!/usr/bin/env python3

import os
import sys
import re
import pathlib
import pyfastx
import pysam
import pandas as pd
import numpy as np
import argparse
from xopen import xopen
from multiprocessing import Pool
from hiseq.fragsize.fragsize import BamFragSizeR1
from hiseq.utils.fastx import Fastx
from hiseq.utils.utils import log, update_obj, run_shell_cmd
from hiseq.utils.file import fx_name, check_path, file_exists, file_prefix
from hiseq.utils.bam import Bam
import hiseq


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
                

def calFreq(x):
    """Calculate the frequency of list
    return dataframe
    index count
    """
    if isinstance(x, list):
        var, freq = np.unique(x, return_counts=True)
        df = pd.DataFrame(data=freq, index=var, columns=['count']).reset_index()
        df.columns = ['length', 'count']
    else:
        df = pd.DataFrame(columns=['length', 'count'])
    return df
    

class BamFragSizeR1(object):
    """Calculate the read size of SE, single BAM

    Parameters
    ----------
    bam: str
        bam file, single

    labels: str
        label

    as_se: bool
        Treat BAM as SE reads

    max_count: int
        maximum number of reads to process, default: [0], all

    strandness: bool
        Check the strandness, default: [False]

    csv_file: str
        File to save the results, in csv format

    sample size = 1000 (SE or PE)

    > Table.csv
    length strand count
    """
    def __init__(self, **kwargs):
        self = update_obj(self, kwargs, force=True)
        self.init_args()


    def init_args(self):
        args_init = {
            'bam': None,
            'outdir': None,
            'labels': None,
            'as_se': False,
            'max_count': 0,
            'strandness': False,
            'csv_file': None,
            'plot_file': None,
        }
        self = update_obj(self, args_init, force=False)
        # only for single bam file
        if not isinstance(self.bam, str):
            raise ValueError('str expected, {} got'.format(
                type(self.bam).__name__))
        if not self.is_bam(self.bam):
            raise Exception('bam={}, is not bam file'.format(self.bam))
        # bai
        if not file_exists(self.bam+'.bai'):
            Bam(self.bam).index()
        # outdir
        if not isinstance(self.outdir, str):
            self.outdir = str(pathlib.Path.cwd())
        check_path(self.outdir)
        # check strandness
        # as_se: paired end only
        if not (self.is_paired(self.bam) and not self.as_se):
            self.as_se = True
        if not self.labels:
            self.labels = file_prefix(self.bam)
        # output files
        if self.csv_file is None:
            self.csv_file = os.path.join(self.outdir, self.labels + '.fragsize.csv')
        if self.plot_file is None:
            self.plot_file = os.path.join(self.outdir, self.labels + '.fragsize.pdf')


    def is_empty(self, bam):
        """Check input file is empty or not
        pysam.Samfile().count
        pysam.Samfile().mapped
        """
        pysam.index(bam)
        sam = pysam.Samfile(bam)
        return sam.count() == 0 and sam.mapped == 0


    def is_bam(self, bam):
        """Check input is BAM file
        update: self.bam

        string
        *.bam
        """
        if bam is None:
            bam = self.bam
        out = False
        if isinstance(self.bam, str):
            out = self.bam.endswith('.bam') and os.path.exists(self.bam)
        return out


    def is_paired(self, bam, topn=1000):
        """Check input bam is Paired end alignment"""
        out = False
        if self.is_bam(bam):
            samfile = pysam.AlignmentFile(bam)
            out = all([read.is_paired for read in samfile.head(topn)])
            samfile.close()
        return out


    def cal_freq(self, x, return_dataframe=True):
        """Calculate the frequency of list
        ['length', 'strand']
        return dataframe
        """
        header = ['length', 'strand', 'count']
        if isinstance(x, list):
            df = pd.DataFrame(x, columns=header).groupby(['length', 'strand']).count().reset_index()
        else:
            df = pd.DataFrame(columns=header)
        if not self.strandness:
            df['strand'] = '*'
        return df


    def cal_frag_size(self, bam=None, chunk=1000000):
        """Extract the read length
        length count id
        """
        if bam is None:
            bam = self.bam
        # for SE or PE
        if self.is_paired(bam):
            pass
        else:
            pass
        # empty check
        if self.is_empty(bam):
            log.error('bam is empty: {}'.format(bam))
            return pd.DataFrame(columns=['length','strand','count','id']) #!!
        # read sam/bam file
        sam = pysam.AlignmentFile(bam)
        counter  = 0
        frag_size_list = []
        frames = []
        for read in sam:
            n_count = parse_fx_id(read.query_name) # !!!
            if self.as_se:
                # reads sizes
                if not read.is_unmapped \
                and not read.is_duplicate > 0:
                    counter += 1
                    strand = '-' if read.is_reverse else '+'
                    frag_size_list.append([read.infer_query_length(), strand, n_count])
            else:
                strand = '*'
                # fragment sizes
                if read.is_proper_pair \
                and not read.is_unmapped \
                and not read.mate_is_unmapped \
                and not read.is_read1 \
                and not read.is_duplicate \
                and read.template_length > 0:
                    counter += 1
                    frag_size_list.append([read.template_length, strand, n_count])
            # sample size
            if self.max_count > 0 and counter  > self.max_count:
                log.info('stop at: {}'.format(counter))
                break # stop
            # chunk
            if counter > 0 and counter%chunk == 0:
                frames.append(self.cal_freq(frag_size_list))
                frag_size_list = [] # empty
                log.info('{} : {} {}'.format('Processed', counter , self.labels))
        # last chunk
        if len(frag_size_list) > 0:
            frames.append(self.cal_freq(frag_size_list))
            frag_size_list = [] # empty
            log.info('{} : {} {}'.format('Processed', counter , self.labels))
        # overall
        df = pd.concat(frames, axis=0).groupby(['length', 'strand']).sum().reset_index()
        df['id'] = self.labels
        return df


    def distribution(self):
        """Basic statistics values
        value + freq
        mean, medium, mode, std, min, max, Q1, Q2, Q3
        """
        if self.freq_table.shape[0] == 0:
            out = pd.DataFrame(
                columns=['mean', 'median', 'mode', 'std', 'min', 'max', 'Q1',
                         'Q2', 'Q3'])
        else:
            val = self.freq_table['length']
            freq = self.freq_table['count']
            inserts = np.repeat(val, freq)
            # statistics
            q_mean = np.mean(inserts)
            q_median = np.median(inserts)
            q_median_dev = np.median(np.absolute(inserts - q_median))
            q_mode = val[np.argmax(freq)]
            q_std = np.std(inserts)
            q_min = np.min(inserts)
            q_max = np.max(inserts)
            q_qual = np.quantile(inserts, [0.25, 0.5, 0.75], axis=0)
            # core distribution
            s = np.array([q_mean, q_median, q_mode, q_std, q_min, q_max]).round(2)
            s = np.append(s, q_qual)
            # DataFrame
            out = pd.DataFrame(s).T
            out.columns = ['mean', 'median', 'mode', 'std', 'min', 'max', 'Q1',
                           'Q2', 'Q3']
        return out


    def _tmp(self):
        """
        Create a tmp file to save json object
        """
        tmp = tempfile.NamedTemporaryFile(prefix='tmp', suffix='.csv',
            delete=False)
        return tmp.name


    def save_as(self, csv_file=None):
        """Save to file"""
        if csv_file is None:
            csv_file = self.csv_file # default
        if csv_file is None:
            csv_file = self._tmp()
        log.info('saving to file: {}'.format(csv_file))
        try:
            self.freq_table.to_csv(csv_file, index=False)
            # save statistics
            stat_file = csv_file + '.stat'
            da = self.distribution()
            da.to_csv(stat_file, sep='\t', index=False)
        except:
            log.warning('failed saving file: {}'.format(csv_file))


    def plot(self, plot_file=None):
        """Generate freq table plot
        line plot, hist
        """
        if plot_file is None:
            plot_file = self.plot_file
        hiseq_dir = os.path.dirname(hiseq.__file__)
        frag_plot_r = os.path.join(hiseq_dir, 'bin', 'qc_fragsize.R')
        stdout = os.path.join(self.outdir, self.labels+'.fragsize.plot.stdout')
        stderr = os.path.join(self.outdir, self.labels+'.fragsize.plot.stderr')
        cmd_file = os.path.join(self.outdir, self.labels+'.fragsize.plot.cmd.sh')
        cmd = ' '.join([
            'Rscript',
            frag_plot_r,
            plot_file,
            self.csv_file,
            '1> {}'.format(stdout),
            '2> {}'.format(stderr)
        ])
        with open(cmd_file, 'wt') as w:
            w.write(cmd+'\n')
        try:
            run_shell_cmd(cmd)
        except:
            log.warning('fragsize.py failed')
        return(plot_file)


    def run(self):
        self.freq_table = self.cal_frag_size(self.bam) # dataframe
        self.save_as() # to csv
        self.plot() # to pdf


def get_args():
    parser = argparse.ArgumentParser(description='hiseq fragsize.py [-o out.csv] <in.fq>')
    parser.add_argument('-o', '--csv-file', dest='csv_file', help='csv file')
    parser.add_argument('bam', help='bam file')
    return parser


def main():
    args = get_args().parse_args()
    args_local = {
        'bam': args.bam,
        'outdir': fx_name(args.bam),
        'csv_file': args.csv_file,
        'as_se': True,
        'strandness': True,
    }
    BamFragSizeR1(**args_local).run()

if __name__ == '__main__':
    main()

#