
# Getting Started


## 1. Run piRNA analysis in one step  

```
# for single file
$ python piRNA_pipe.py -i demo.fq.gz -o output --trimmed -g dm6 -p 8

# for multiple files
$ python piRNA_pipe.py -i demo1.fq.gz demo2.fq.gz -o output --trimmed -g dm6 -p 8 -j 2
```

output file structure: 

```
report/
├── 00.raw_data
├── 01.clean_data
├── 02.collapse
├── 03.overlap
├── 04.smRNA
├── 05.miRNA
├── 06.size_exclude
├── 06.size_select
├── 07.te
├── 08.piRNA_cluster
├── 09.genome
├── 10.unmap
├── 11.stat
├── 12.report
├── 13.genome
└── config
```

`11.stat/fx_stat.seqs.csv` - statistics of sequences  
`11.stat/fx_stat.reads.csv` - statistics of reads  
`12.report/smRNA_report.html` - final report in HTML format  


## 2. count reads for fastq file

**Important**: fastq files were collapsed after `02.collapse` step, in this pipeline. the read name were in `{num1}-{num2}` format:

`num1` - indicate the ranking of the sequence in the file  
`num2` - indicate the number of the sequence (read count)

Using the script `fx_count.py` to count reads for fastq file: 

```
$ python fx_count.py -o fx_stat.csv demo.fq.gz

demo    839312  2517074
```

`demo` - name of the fastq file 
`839312` - number of sequences (unique sequence) 
`2517074` - number of reads 


## 3. check fragment size for fastq file 

see above section, need to parse the read count from read name.

```
$ python fx_fragsize.py -o out.csv demo.fq.gz 
```


## 4. check fragment size for bam file 

see section2, parse the read count from read name 


```
$ python bam_fragsize.py -o out.csv demo.bam
```

## 5. check base-content for fastq file

see section2, parse the read count from read name 

```
$ 
```




