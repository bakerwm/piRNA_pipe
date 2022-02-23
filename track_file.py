import os
import sys
import pathlib
import numbers
import pysam
import pyBigWig
from hiseq.utils.helper import *
from utils import *

"""Processing the track files: bigWig, BED, ...
Using the following packages:

CLI tool: [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks)
API tool: [pyBigWig](https://github.com/deeptools/pyBigWig), 
          [pybedtools](https://github.com/daler/pybedtools)
          
Example:
# for all chromosomes
>>> get_tracks(bw_list, outdir=out)

# for BED regions
>>> get_tracks(bw_list, outdir=out, region_list=in_bed)

# specify colors
>>> get_tracks(bw_list, outdir=out, region_list=in_bed, colors=['red', 'green', 'blue'])

"""

################################################################################
## main functions bigwig tracks
################################################################################    
def bw_to_config(x, **kwargs):
    """Generate the config.ini for single file

    Parameters
    ----------
    x : str
        path to the bigwig file

    Keyword arguments
    -----------------
    color: color for the file, default '#333333' grey
    min_value: display on y-axis, default: 0
    max_value: display on y-axis, default: 'auto'
    height: the height for the track, defautl: 2
    nbins: number of bins, default: 700

    Example output
    --------------
    [track-01]
    file = 
    title = 
    height = 2
    color = 
    min_value = auto
    max_value = auto
    number_of_bins = 700
    file_type = bigwig
    """
    if not isinstance(x, str) or not file_exists(x):
        log.error('x={} not valid or not exists'.format(x))
        return None
    # default
    color = kwargs.get('color', '#333333')
    y_min = kwargs.get('min_value', 0) 
    y_max = kwargs.get('max_value', 'auto')
    height = kwargs.get('height', 2)
    nbins = kwargs.get('nbins', 700)
    xname, xtype = os.path.splitext(os.path.basename(x))
    xtype = xtype.lstrip('.')
    return '\n'.join([
        '[track-{}]'.format(xname),
        'file = {}'.format(x),
        'title = {}'.format(xname),
        'height = {}'.format(height),
        'color = {}'.format(color),
        'min_value = {}'.format(y_min),
        'max_value = {}'.format(y_max),
        'number_of_bins = {}'.format(nbins),
        'file_type = {}'.format(xtype.lower()),
    ])  

        
def get_track(x, chr, **kwargs):
    """Generate the config.ini and run.sh command for single region
    
    Parameters
    ----------
    x : list
        list of bigWig files

    chr : str
        specify the chr name

    Keyword arguments
    -----------------    
    outdir: Directory to save the config.ini and plot
            if `None`, set the `cwd()`.
    start:  Starting position
    end:    Ending position
    y_min:  The minimum value on y-axis, default: 0
    y_max:  The maximum value on y-axis, default: 'auto'
    colors: Set colors for the tracks/bigWig files, 
            if `auto`, use colors by `fishualize` package
            from R, or set a list of colors by name, 
            `['red', 'green', 'blue']`, the number should
            match the bigWig files
    config: The name of `ini` file, defualt: `outdir/config.ini`
    plot_file: The name of the plot file, 'chr_start_end.png'
        
    Example
    [track-01]
    ...
    [spacer]
    [x-axis]
    """
    ###########################################################################
    ## Config check
    ## supported files
    if isinstance(x, list):
        x = [b for b in x if file_exists(b)] # remove None
    else:
        log.error('x expect list, got {}'.format(type(x).__name__))
        return None
    if len(x) == 0:
        log.error('x, no files valid')
        return None
    x = [file_abspath(i) for i in x] # convert to absolute_path
    ## arguments
    chr_list = load_regions(x[0]) # dict: {chr:(start, end)}
    if not chr in chr_list:
        log.error('chr={} no valid, for example: {}'.format(chr, chr_list[:5]))
        return None                  
    outdir = kwargs.get('outdir', None)
    if not isinstance(outdir, str):
        outdir = str(pathlib.Path.cwd())
    outdir = file_abspath(outdir)
    project_dir = os.path.join(outdir, chr)
    check_path(project_dir)
    ## start/end
    start = kwargs.get('start', 1)
    end = kwargs.get('end', 1000)
    ## color
    i_colors = kwargs.get('colors', None)
    if isinstance(i_colors, str): 
        x_colors = i_colors * len(x)
    elif isinstance(i_colors, list):
        if len(i_colors) >= len(x):
            x_colors = i_colors[:len(x)]
        else:
            x_colors = fish_colors(len(x))
    else:
        x_colors = fish_colors(len(x))
    ## y-axis
    y_max = kwargs.get('y_max', 'auto')
    if y_max == 'auto' or isinstance(y_max, numbers.Number):
        pass
    else:
        y_max = get_bw_y_axis(x, chr, start=start, end=end)[1] # (y_min, y_max)
    y_min = kwargs.get('y_min', 0)
    if y_min == 'auto' or isinstance(y_min, numbers.Number):
        pass
    else:
        y_min = 0
    ###########################################################################
    ## generate ini
    args = kwargs.copy()
    args.update({
        'start': start,
        'end': end,
        'min_value': y_min,
        'max_value': y_max,
    })
    msg_list = []
    for bw,color in zip(x, x_colors):
        args.update({'color': color})
        msg = bw_to_config(bw, **args)
        if msg:
            msg_list.append(msg)
    ## add extra
    msg_list += ['[spacer]', '[x-axis]']
    config_text = '\n'.join(msg_list)
    ## save to file
    config_name = kwargs.get('config', 'config.ini')
    config_file = os.path.join(project_dir, config_name)
    with open(config_file, 'wt') as w:
        w.write(config_text+'\n')
    ###########################################################################
    ## generate command-line
    plot_stdout = os.path.join(project_dir, 'stdout.log')
    plot_stderr = os.path.join(project_dir, 'stderr.log')
    plot_cmd = os.path.join(project_dir, 'run.sh')
    plot_file = kwargs.get('plot_file', None)
    if isinstance(plot_file, str):
        check_path(os.path.dirname(plot_file))
        pass
    else:
        plot_file = os.path.join(project_dir, 
                                 '{}_{}_{}.png'.format(chr, start, end))
    ## command
    cmd = ' '.join([
        '{}'.format(shutil.which('pyGenomeTracks')),
        '--tracks {}'.format(config_file),
        '--region {}:{}-{}'.format(chr, start, end),
        '-o {}'.format(plot_file),
        '--dpi 150',
        '--trackLabelFraction 0.2',
        '1> {}'.format(plot_stdout),
        '2> {}'.format(plot_stderr),
    ])
    with open(plot_cmd, 'wt') as w:
        w.write(cmd+'\n')
    if file_exists(plot_file):
        log.info('get_track() skipped, file exists: {}'.format(plot_file))
    else:
        os.system(cmd)
    return plot_file



def get_tracks(x, **kwargs):
    """Generate the plots

    Parameters
    ----------
    x : list
        list of bigWig files
        
    Keyword parameters
    ------------------
    outdir: Directory to save the config.ini and plot
            if `None`, set the `cwd()`.
    region_list: regions in BED format
    y_min:  The minimum value on y-axis, default: 0
    y_max:  The maximum value on y-axis, default: 'auto'
    colors: Set colors for the tracks/bigWig files, 
            if `auto`, use colors by `fishualize` package
            from R, or set a list of colors by name, 
            `['red', 'green', 'blue']`, the number should
            match the bigWig files
    nbins:  Number of bins, default: 700
    height: The height of each track
    config: The name of `ini` file, defualt: `outdir/config.ini`
    plot_frefix: The name of the plot file, 'chr_start_end.png'
    
    >>> get_config
    
    """
    if isinstance(x, list):
        x = [b for b in x if file_exists(b)] # remove None
    else:
        log.error('x expect list, got {}'.format(type(x).__name__))
        return None
    if len(x) == 0:
        log.error('x, no files valid')
        return None
    x = [file_abspath(i) for i in x] # convert to absolute_path
    region_list = kwargs.get('region_list', None)
    r = region_list if file_exists(region_list) else x[0]
    r_list = load_regions(r)
    ## config
    msg = '\n'.join([
        '-'*80,
        'track_file: {}'.format(',\n'.join(x)),
        'outdir: {}'.format(kwargs.get('outdir', '')),
        'region_list: {}'.format(r),
        '-'*80,
    ])
    print(msg)    
    ## run
    i = 0
    ia = len(r_list)
    for chr,pos in r_list.items():
        i += 1
        args_k = {
            'start': pos[0],
            'end': pos[1],
        }
        kwargs.update(args_k)
        print('{}/{} - {}'.format(i, ia, chr))
        p = get_track(x, chr, **kwargs)        
        # break


################################################################################
## Preparing the track files: bigWig, ...
################################################################################
def bam_to_bw(bam, bw, scale=1, gsize=0, threads=8):
    """Convert bam to bigWig, using bamCoverage (deeptools)
    Parameters
    ----------
    bam : str
        Path to the bam file
        
    bw : str
        Path to the bigWig file
        
    scale : float
        The scale factor, default: 1
        
    gsize : int
        The genome size, if `0`, extract from the BAM header.
        
    >>> bam_to_bw(in_bam, out_bw, 1)
        
    """
    # check bam exists
    if isinstance(bam, str) and file_exists(bam):
        pass
    else:
        log.error('bam= not exists, {}'.format(bam))
        return None
    # bw directory
    if not isinstance(bw, str):
        log.error('bw= expect str, got {}'.format(type(bw).__name__))
        return None
    bw_dir = os.path.dirname(bw)
    check_path(bw_dir)
    # get the gsize
    ref_size = gsize if gsize > 0 else sum(pysam.AlignmentFile(bam).header.lengths)
    # cmd
    cmd_sh = os.path.join(os.path.dirname(bam), os.path.basename(bam)+'.bam2bw.sh')
    cmd_log = os.path.join(os.path.dirname(bam), os.path.basename(bam)+'.bam2bw.log')
    cmd = ' '.join([
        '{}'.format(shutil.which('bamCoverage')),
        '-b {}'.format(bam),
        '-o {}'.format(bw),
        '--scaleFactor {}'.format(scale),
        '--effectiveGenomeSize {}'.format(ref_size),
        '-p {}'.format(threads),
        '--binSize 1 --normalizeUsing CPM',
        '2> {}'.format(cmd_log),
    ])
    with open(cmd_sh, 'wt') as w:
        w.write(cmd+'\n')
    # run
    if os.path.exists(bw):
        print('bam_to_bw() skipped, file exists: {}'.format(bw))
    else:
        os.system(cmd)
    return bw
    


    
def run_bam_to_bw(x, group='te', strandness=None, unique='unique', outdir=None):
    """Convert bam to bigwig for specific group, in pipe() project
    
    Parameters
    ----------
    x : str
        path to the project_dir (pipe)

    group : str
        The group of reads, choose from ['smRNA', 'miRNA', 'te', 'piRC', 'genome'], 
        default: `te`

    strandness : str
        Extract reads from the strand: ['fwd', 'rev', 'both', None], default: 'both'        

    unique : str
        The unique or multi mapping, ['unique', 'multi', 'both'], default: 'unique'
        
    outdir : str
        The directory saving the bigWig files, default `None`, the same directory 
        as the bam file.
    """
    if not is_pipe_dir(x):
        log.error('x is not pipe() project dir, {}'.format(x))
        return None
    group_list = ['smRNA', 'miRNA', 'te', 'piRC', 'genome', 'genome2']
    if not group in group_list:
        log.error('group={} not valid, choose: {}'.format(group, group_list))
        return None
    strand_list = ['fwd', 'rev', 'both', None, '1', '-1', '0', '*']
    if not strandness in strand_list:
        log.error('strandness={} not valid, choose: {}'.format(strandness, strand_list))
        return None
    unique_list = ['unique', 'multi', 'both']
    if not unique in unique_list:
        log.error('unqiue={} not valid, choose: {}'.format(unique, unique_list))
    # get the bam file
    bam = get_x_file(x, 'bam', group, unique)
    bw = get_x_file(x, 'bigwig', group, unique, check_exists=False)
    print('!AAAA-2', bam, bw)
    prefix = os.path.splitext(bam)[0]
    if not file_exists(bam): # exists
        log.error('bam file not exists: {}'.format(bam))
        return None
    # scale
    n_map = get_x_map(x, 'map') # clean - unmap
    # scale = 1e6/n_map if n_map > 0 else 1
    scale = 1
    if strandness in ['fwd', '1', 'both']:
        # forward
        bam_fwd = prefix+'.fwd.bam'
        bw_fwd = prefix+'.fwd.bigwig'
        if not file_exists(bam_fwd):
            pysam.view('-bhS', '-F', '16', '-o', bam_fwd, bam, catch_stdout=False)
            pysam.index(bam_fwd)
        bam_to_bw(bam_fwd, bw_fwd, scale=scale)
    if strandness in ['rev', '-1', 'both']:
        # reverse
        bam_rev = prefix+'.rev.bam'
        bw_rev = prefix+'.rev.bigwig'
        if not file_exists(bam_rev):
            pysam.view('-bhS', '-f', '16', '-o', bam_rev, bam, catch_stdout=False)
            pysam.index(bam_rev)
        bam_to_bw(bam_rev, bw_rev, scale=scale)
        out = [bw_fwd, bw_rev]
    if strandness in ['0', '*', None]:
        bam_to_bw(bam, bw, scale=scale)
    # ouptut
    if strandness in ['both']:
        out = [bw_fwd, bw_rev]
    elif strandness in ['fwd', '1']:
        out = [bw_fwd, None]
    elif strandness in ['rev', '-1']:
        out = [bw_rev, None]
    else:
        out = [bw, None]
    return out




################################################################################
## functions for pyGenomeTracks
################################################################################
def load_regions(x):
    """Extract the regions from file
    support: bam/bigwig/bed
    
    Parameters
    ----------
    x : str
        The input file, in ['bam', 'bigwig', 'bed'] format, 
        if `bam` or `bigwig`, return the whole chromosome as region
        or specify the regions in BED format, (BED, -1 start)
        
    output format:
    {chr:(start, end)}
    """
    if not isinstance(x, str) or not file_exists(x):
        log.error('x={} expect str, got {}'.format(type(x).__name__))
        return None
    x_type = os.path.splitext(x)[1]
    x_type = x_type.lower()
    ft_list = ['.bam', '.bw', '.bigwig', '.bed', '.narrowpeak']
    if not x_type in ft_list:
        log.error('x extsion: {} not valid, choose: {}'.format(x_type, ft_list))
    ## for each type of file
    out = None
    if x_type == '.bam':
        try:
            h = pysam.AlignmentFile(x).header
            out = {i:(1, k) for i,k in zip(h.references, h.lengths)}
        except:
            log.error('reading file failed, {}'.format(x))
    elif x_type in ['.bw', '.bigwig']:
        try:
            bw = pyBigWig.open(x)
            d = bw.chroms()
            out = {k:(1, v) for k,v in d.items()}
            bw.close()
        except:
            log.error('reading file failed, {}'.format(x))
    elif x_type in ['.bed', '.bed3', '.bed6', '.bed12', '.narrowpeak']:
        try:
            out = {}
            for bed in pybedtools.BedTool(x):
                out[bed.chrom] = (int(bed.start)+1, int(bed.end))
        except:
            log.error('reading file failed, {}'.format(x))
    return out
            

def fish_colors(n=4, colorPal=1):
    """Pick colors

    Parameters
    ----------
    n : int
        number of colors return, default: 4
        
    colorPal : int or str
        choose the group of colors, 1 to N, or the name of fish
        candidate list: [1, 2, 3, 'Scarus_hoefleri', 'Gramma_loreto',
        'Centropyge_loricula'], default: [1]

    # colors from fishualize R package
    # https://nschiett.github.io/fishualize/articles/overview_colors.html

    Colors from R package: 
    > fishualize::fish_pal(option = "Scarus_hoefleri")(10)
     [1] "#D2372CFF" "#E25719FF" "#F0780BFF" "#F89814FF"
     [5] "#EFB52BFF" "#C6CB43FF" "#7AD45CFF" "#00CE7DFF"
     [9] "#00BCABFF" "#0499EAFF"

    > fishualize::fish_pal(option = "Gramma_loreto")(10)
     [1] "#020122FF" "#1E2085FF" "#4029CBFF"
     [4] "#6628EEFF" "#901CEDFF" "#B804CAFF"
     [7] "#D61693FF" "#E6445DFF" "#EE7A30FF"
    [10] "#F0BF0BFF"

    > fishualize::fish_pal(option = "Centropyge_loricula")(10)
     [1] "#8F1D1EFF" "#B30029FF" "#DF002AFF"
     [4] "#FF7D1AFF" "#FFBD17FF" "#E7BE5AFF"
     [7] "#988591FF" "#0043A0FF" "#001A72FF"
    [10] "#000000FF"
    """
    # pre-defined colors
    c1 = ['#D2372C', '#E25719', '#F0780B', '#F89814', '#EFB52B', 
        '#C6CB43', '#7AD45C', '#00CE7D', '#00BCAB', '#0499EA']
    c2 = ['#020122', '#1E2085', '#4029CB', '#6628EE', '#901CED', 
        '#B804CA', '#D61693', '#E6445D', '#EE7A30', '#F0BF0B']
    c3 = ['#8F1D1E', '#B30029', '#DF002A', '#FF7D1A', '#FFBD17', 
        '#E7BE5A', '#988591', '#0043A0', '#001A72', '#000000']
    # RGB (10 colors, repeat twice, 20 items)
    color_d = {
        'Scarus_hoefleri': c1*2,
        'Gramma_loreto': c2*2,
        'Centropyge_loricula': c3*2,
    }
    # get the fish_list
    fish_list = list(color_d.keys())
    # determine the fish (colorBy)
    if isinstance(colorPal, int):
        colorPal = colorPal - 1 # to 0-indexed
        if not colorPal in range(len(fish_list)):
            colorPal = 0
        fish = fish_list[colorPal]
    elif isinstance(colorPal, str):
        fish = colorPal if colorPal in color_d else 'Scarus_hoefleri'
    else:
        fish = 'Scarus_hoefleri'
    # output colors
    return color_d.get(fish, c1)[:n]


def bw_stats(x, chr=None, **kwargs):
    """Get the stats for bigWig file, in specific region
    
    Parameters
    ----------
    x : str
        bigwig file
        
    chr : str 
        The chromosome name, if `None`, pick the first chromosome

    Keyword arguments
    -----------------
    start: Starting position
    end:   Ending position
    type:  Summary type (mean, min, max, coverage, std), default 'mean'.
    nBins: Number of bins into which the range should be divided before
           computing summary statistics. The default is 1.
    exact: By default, pyBigWig uses the same method as Kent's tools from UCSC
           for computing statistics. This means that 'zoom levels' may be
           used, rather than actual values (please see the pyBigWig repository
           on github for further information on this). To avoid this behaviour,
           simply specify 'exact=True'. Note that values returned will then
           differ from what UCSC, IGV, and similar other tools will report.

    >>> import pyBigWig
    >>> bw = pyBigWig.open("test/test.bw")
    >>> bw.stats("1", 0, 3)
    [0.2000000054637591]
    pyBigWig.open(bw).stats('Idefix', 1, 7411, type='max')    
    """
    bw = pyBigWig.open(x)
    if chr is None: # choose chromosome
        chr = list(bw.chroms().keys()).pop(0)
    s = bw.stats(chr, **kwargs)
    bw.close()
    return s


def get_bw_y_axis(x, chr=None, start=0, end=0, fix_end=True):
    """Determine the y-min, y-max, for bigwig files in specific region
    
    Parameters
    ----------
    x : list
        A list of bigwig files
        
    chr : str
        The name of chromosome, if `None`, choose the first chr 
        
    start : int 
        The start position 
        
    end : int
        The end position, 0 indicate the end of the chromosome
        
    fix_end : bool
        Fix the value, to nearest int (by breaks: 5)
        
    Output
    ------
    The format for pyGenomeTracks ini file, float or 'auto'
    out = (ymin, ymax)
    """
    if not isinstance(x, list):
        log.error('x, expect list, got {}'.format(type(x).__name__))
        return ('auto', 'auto')
    x = [i for i in x if file_exists(i)]
    if len(x) == 0:
        log.error('x, no files')
        return ('auto', 'auto')
    # fix start, end
    args = {}
    if start > 0:
        args['start'] = start
    if end > 0:
        args['end'] = end
    ## min
    args_min = args.copy()
    args_min['type': 'min']
    s = [bw_stats(bw, chr=chr, **args_min) for bw in x]
    smin = round(min(s).pop(0), 1)
    if fix_end:
        breaks = 5 # 3, 5, 10, ...
        smin = (int(smin/breaks))*breaks
    ## max
    args_max = args.copy()
    args_max['type': 'max']
    s_max = [bw_stats(bw, chr=chr, **args_min) for bw in x]
    smax = round(max(s).pop(0), 1)
    if fix_end:
        breaks = 5 # 3, 5, 10, ...
        smax = (int(round(smax/breaks, 1))+1)*breaks
    return (smin, smax)
    


    
    
    

def main():
    if len(sys.argv) < 3:
        print('Usage: track_file.py <prj.list> <outdir>')
        sys.exit(1)
    # args
    prj_txt = sys.argv[1]
    outdir = sys.argv[2]
    with open(prj_txt) as r:
        prj_list = r.readlines()
    prj_list = [i.strip() for i in prj_list]
    bw_list = [get_x_file(i, 'bigwig', 'te') for i in prj_list]
    get_tracks(bw_list, outdir=outdir)
#     # bw_to_tracks(bw_list, outdir=outdir)
    
    
if __name__ == '__main__':
    main()
    
    
