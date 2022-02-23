#!/bin/bash
# calculate base-content for fastq or fasta files

# Author: Ming Wang
# Date: 2018-05-22

# functions
function is_gzip() {
    if [[ -z $(gzip -t $1 2>&1) ]] ; then
        echo "yes"
    else
        echo "no"
    fi
}

function seq_type() {
    fn=$1
    [[ $(is_gzip $1) = "yes" ]] && fn_viewer="zcat" || fn_viewer="cat"
    stat=$(${fn_viewer} ${fn} | \
        head -n 1000 | \
        awk '{if(NR%4==1){t=substr($1,1,1); print t}}' | \
        uniq)
    #echo ${stat}
    case ${stat} in 
        \@)    seq_tag='fastq' ;;
        \>)    seq_tag='fasta' ;;
        *)    echo2 "unknown sequence type" "error"
    esac
    echo ${seq_tag}
}

function awk_fastq() {
    # save to two files
    # len_dis
    # base_content
    cat <<EOF
{
    if(NR%4==2) {
        lenmax=length(\$1)
        len[lenmax]++
        if (max == 0) { max=lenmax }
        for (i=1; i<=max; i++) {
            q=substr(\$1,i,1)
            s[i][q]++
            #print c
        }
    }
} 
END {
    print "position","A","C","G","T","N" > f1
    for (i in s) {
        if (isarray(s[i])) {
            if (s[i]["A"] == "") {s[i]["A"]=0}
            if (s[i]["C"] == "") {s[i]["C"]=0}
            if (s[i]["G"] == "") {s[i]["G"]=0}
            if (s[i]["T"] == "") {s[i]["T"]=0}
            if (s[i]["N"] == "") {s[i]["N"]=0}
            print i,s[i]["A"],s[i]["C"],s[i]["G"],s[i]["T"],s[i]["N"] > f1
        }
    }
    for (j in len) {
        print j,len[j] > f2
    }
}
EOF
}

function awk_fasta() {
    # save to two files
    # len_dis
    # base_content
    cat <<EOF
{
    if(NR%2==1) {
        lenmax=length(\$1)
        len[lenmax]++
        if (max == 0) { max=lenmax }
        for (i=1; i<=max; i++) {
            q=substr(\$1,i,1)
            s[i][q]++
            #print c
        }
    }
} 
END {
    print "position","A","C","G","T","N" > f1
    for (i in s) {
        if (isarray(s[i])) {
            if (s[i]["A"] == "") {s[i]["A"]=0}
            if (s[i]["C"] == "") {s[i]["C"]=0}
            if (s[i]["G"] == "") {s[i]["G"]=0}
            if (s[i]["T"] == "") {s[i]["T"]=0}
            if (s[i]["N"] == "") {s[i]["N"]=0}
            print i,s[i]["A"],s[i]["C"],s[i]["G"],s[i]["T"],s[i]["N"] > f1
        }
    }
    for (j in len) {
        print j,len[j] > f2
    }
}
EOF
}

function content_fastx() {
    # calculate the base content for fastx file
    fn=$1 # fastx file
    f1=$2 # stat_txt
    f2=$3 # lendis_txt
    [[ -z $2 ]] && echo "require stat.txt file" && exit 1
    [[ -z $3 ]] && echo "require lendis.txt file" && exit 1
    [[ -z $4 ]] && n=0 || n=$4
    fs=$(tempfile -s ".awk")
    # >&2 echo $fs
    if [[ $(seq_type $fn) = "fasta" ]] ; then
        awk_fasta > $fs
    elif [[ $(seq_type $fn) = "fastq" ]] ; then
        awk_fastq > $fs
    else
        echo "unknown file: $fn"
        exit 1
    fi
    ##
    [[ $(is_gzip $fn) = "yes" ]] && fn_viewer="zcat" || fn_viewer="cat"
    ${fn_viewer} $fn | awk -v max=$n -v f1=$f1 -v f2=$f2 -f $fs
    [[ -f $fs ]] && rm $fs
}

function get_rscript() {
    cat <<EOF
#!/usr/bin/env Rscript
## create base-content figure
## Author: Ming Wang
## Date: 2018-05-20

args = commandArgs(trailingOnly = TRUE)

# test arguments
if (length(args) != 2) {
  stop("Rscript base_content_figure.R <data.dir> <out.prefix>", call. = FALSE)
}
indir <- args[1] # statdir
outname <- args[2] # name

f1 <- file.path(indir, "base_content.txt")
f2 <- file.path(indir, "length_distribution.txt")

if (! file.exists(f1) & file.exists(f2)) {
  stop("file not exists: base_content.txt, length_distribution.txt")
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

##----------------------------------------------------------------------------##
df <- read.table(f1, header = TRUE, sep = " ")
df[is.na(df)] <- 0 # NA to 0
df <- tidyr::gather(df, base, count, -position)
df <- filter(df, base %in% c("A", "C", "G", "T")) %>%
  mutate(base = factor(base, levels = c("A", "C", "G", "T"))) %>%
  group_by(position) %>%
  mutate(freq = count / sum(count))

p1 <- ggplot(df, aes(position, freq, fill = base)) +
  geom_bar(position = "fill", stat = "identity") +
  xlab("Position in read (bp)") +
  ylab("Per base sequence content (%)") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  scale_fill_manual(values = c("green3", "blue3", "grey30", "red3")) +
  theme_classic()

p2 <- ggplot(df, aes(position, freq, color = base)) +
  geom_line(size = .5) +
  #geom_bar(position = "fill", stat = "identity") +
  xlab("Position in read (bp)") +
  ylab("Per base sequence content (%)") +
  scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 100, 20)) +
  scale_color_manual(values = c("green3", "blue3", "grey30", "red3")) +
  theme_classic()

##----------------------------------------------------------------------------##
df2 <- read.table(f2, header = FALSE, sep = " ")
names(df2) <- c("length", "count")
# if only one record in df2
if (nrow(df2) == 1) {
  n <- df2[1, "length"]
  df2x <- data.frame(length = c(n - 1, n + 1),
                     count = 0)
  df2 <- rbind(df2, df2x)
}

p3 <- ggplot(df2, aes(length, count)) +
  geom_line(size = .5, color = "red3") +
  xlab("Length of reads (nt)") +
  ylab("Number of reads") +
  theme_classic()

##----------------------------------------------------------------------------##
png1 <- file.path(indir, paste0(outname, ".base_content.bar.png"))
png2 <- file.path(indir, paste0(outname, ".base_content.line.png"))
png3 <- file.path(indir, paste0(outname, ".length_distribution.line.png"))
ggsave(png1, p1, width = 5, height = 3, units = "in")
ggsave(png2, p2, width = 5, height = 3, units = "in")
ggsave(png3, p3, width = 5, height = 3, units = "in")
##----------------------------------------------------------------------------##

EOF
}

function make_plot() {
    # make plot using R script
    txt=$1
    png=$2
    [[ -z $1 ]] && echo "require stat.txt file, exiting ..." && exit 1
    [[ -z $2 ]] && echo "require out.png file, exiting ..." && exit 1
 
    # tempfile
    rfile=$(tempfile -s ".R")
    get_rscript > $rfile

    # make plot
    [[ -f $rfile ]] && Rscript $rfile $txt $png

    # remove tempfile
    [[ -f $rfile ]] && rm $rfile
}

## main
[[ $# -lt 2 ]] && echo "Usage: bash $(basename $0) <in.fa|fq> <outdir> <N-bases>" && exit 0
f=$1
o=$2
[[ -z $3 ]] && nbases=0 || nbases=$3
[[ ! -d $o ]] && mkdir -p ${o}
#echo "output directory not exists, exiting ..." && exit 1

stat_txt="${o}/base_content.txt"
lendis_txt="${o}/length_distribution.txt"
# content_fastx $f $nbases
content_fastx $f $stat_txt $lendis_txt $nbases 

# make plot using R
outname="myplot"
make_plot $o $outname


## EOF
