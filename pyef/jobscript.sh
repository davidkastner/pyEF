#$ -S /bin/bash          # Specifies interpreting shell for the job
#$ -N pyEF               # Identifier for the job
#$ -l h_rt=48:00:00      # Specify maximum runtime
#$ -l h_rss=40G          # Specify requested memory
#$ -q (gpus|gpusbig)     # Specify which queue to submit the job to
#$ -l gpus=1             # Always set to 1 no matter how many gpus are used.
#$ -pe smp 2             # Number of parallel processes to run (defines how many GPUs)
# -fin run.py            # Copy input file in the running directory
# -fin geometry.py
# -fin analysis.py 

source ~/.bashrc
conda activate molsimp
export OMP_NUM_THREADS=2 # Number of threads to be spawned (needs to be the same as -smp)

echo "running cage jobs!"
cd /home/manets12/DavidMimichrome 
python run.py > Efield_multiwfn_ExampleTest.log
