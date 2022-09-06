#!~/bin/bash

RepeatModeler-2.0.1/BuildDatabase -name chr_rpdb -engine rmblast chr.fa 
RepeatModeler-2.0.1/RepeatModeler -database chr_rpdb -pa 20 
RepeatMasker/RepeatMasker -e rmblast -lib chr_rpdb-families.fa  -pa 20  chr.fa -dir ./ -gff  -norna

