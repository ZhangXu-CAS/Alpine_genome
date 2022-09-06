#! ~/bin/bash
../EvmUtils/partition_EVM_inputs.pl --genome ../chr.fa\
     --gene_predictions ../gene_pred2.gff3 \
	  --transcript_alignments ../trans.aligment.gff3 \
        --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
	 
../EvmUtils/write_EVM_commands.pl --genome ../chr.fa --weights `pwd`/weights.txt \
      --gene_predictions ../gene_pred2.gff3 \
	  --transcript_alignments ../trans.aligment.gff3 \
      --output_file_name evm.out  --partitions partitions_list.out >  commands.list
../EvmUtils/execute_EVM_commands.pl commands.list | tee run.log
../EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out

../EvmUtils/convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome ../chr.fa
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3
