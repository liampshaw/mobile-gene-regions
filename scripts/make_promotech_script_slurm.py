import argparse

# NOTE: not entirely consistent with the modified script I actually use,
# run_promotech_slurm.sh
# but pretty close (just haven't made consistent)

def get_options():
    parser = argparse.ArgumentParser(description='Creates slurm cluster jobs for running PromoTech on multifasta file.')
    parser.add_argument('--fasta_list', help='List of fasta files', type=str)
    return parser.parse_args()

def main():
	args = get_options()
	print("#!/bin/bash")
	print("module load Anaconda3/2020.11")
	print('eval "$(conda shell.bash hook)"')
	print('conda activate promotech_env')
	print('overall_outdir=promoter_predictions')
	print('outdir=tmp-$SLURM_ARRAY_TASK_ID')
	print('fasta=$(sed "${SLURM_ARRAY_TASK_ID}q;d" '+args.fasta_list+")")
	print('echo $name')
	print('time python /well/shaw/users/amu125/programs/PromoTech/promotech.py -pg -m RF-HOT -f $fasta -o $outdir')
	print('time python /well/shaw/users/amu125/programs/PromoTech/promotech.py -g -t 0.6 -i $outdir -o $outdir')
	print('rm $outdir/*data')
	print('mv $outdir/genome_predictions.csv $overall_outdir/"$accession"-promoter_predictions.tsv')
	print('rm -r $outdir/*')
	

if __name__=="__main__":
    main()
