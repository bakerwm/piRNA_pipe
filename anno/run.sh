bowtie -m 1 -x idx/piRC -p 8 -S --un not_piRC.fq --no-unal input.fq.gz > to_piRC.sam 2> to_piRC.log 
bowtie -k 1 -x idx/te   -p 8 -S --un not_te.fq --no-unal   not_piRC.fq > to_te.sam 2> to_te.log 
bowtie -k 1 -x ~/data/genome/dm6/bowtie_index/dm6 -p 12 -S --un not_dm6.fq --no-unal not_te.fq > to_genome.sam 2> to_genome.log
