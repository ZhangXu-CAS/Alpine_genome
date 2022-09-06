#PBS -N Sau_NextDenovo
#PBS -l nodes=1:ppn=20
#PBS -q fat 

###NextDenovo assemble
PBS_O_WORKDIR=/public/home/WHC_zhangx/Saussurea_genome/assembly_files/NextDenovo
cd $PBS_O_WORKDIR
./nextDenovo run.cfg2.1 
echo "nextDenovo done"

###file of nextDenovo input: run.cfg2.1
#[General]
#job_type = local # local, slurm, sge, pbs, lsf
#job_prefix = Sau_nextDenovo2.1
#task = all # all, correct, assemble
#rewrite = yes # yes/no
#deltmp = yes 
#rerun = 3
#parallel_jobs = 15 # number of tasks used to run in parallel
#input_type = raw
#input_fofn = ./input.fofn
#workdir = ./Saussurea_assemble2.1
#
#[correct_option]
#read_cutoff = 1k
#seed_cutoff = 13138 # minimum seed length
#blocksize = 5g
#pa_correction = 15 # number of corrected tasks used to run in parallel, overwrite ${parallel_jobs} only for this step
#seed_cutfiles = 20
#sort_options = -m 50g -t 40 -k 40 # -k, max depth of each overlap, should <= average sequencing depth
#minimap2_options_raw = -x ava-ont -t 40  # change to ava-pb for PacBio CLR data
#correction_options = -p 15
#
#[assemble_option]
#minimap2_options_cns = -x ava-ont -t 40 -k17 -w17 # change to ava-pb for PacBio CLR data
#nextgraph_options = -a 1
