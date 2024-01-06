#$ -S /bin/bash             # Specifies interpreting shell for the job
#$ -N pyEF                  # Identifier for the job
#$ -l h_rt=48:00:00         # Specify maximum runtime
#$ -l h_rss=40G             # Specify requested memory
#$ -q gpusnew               # Specify which queue to submit the job to
#$ -l gpus=1                # Always set to 1 no matter how many gpus are used.
#$ -pe smp 2                # Number of parallel processes to run (defines how many GPUs)
# -fin job_paths.in
# -fin metal_indices.in

export OMP_NUM_THREADS=2    # Number of threads to be spawned (needs to be the same as -smp)
source ~/.bashrc   
conda activate molsimp   

# --geom                      Whether to perform a geometry check
# --esp                       Whether to perform ESP analysis
# --jobs_file                 Paths to all the QM jobs
# --metals_file               The indices of the metals for each job

RUN="/home/kastner/packages/pyEF/pyef/run.py"
JOBS="/home/kastner/packages/pyEF/demo/jobs.in"
METALS="/home/kastner/packages/pyEF/demo/metals.in"
python $RUN --geom --esp --jobs_file $JOBS --metals_file $METALS  > pyEF.log

