#$ -S /bin/bash      # specifies interpreting shell for the job
#$ -N MinimiChrome_Analysis            # identifier for the job
#$ -l h_rt=48:00:00 # specify maximum runtime
#$ -l h_rss=40G       # specify requested memory
#$ -q gpus|gpusbig           # specify which queue to submit the job to
#$ -l gpus=1         # Always set to 1 no matter how many gpus are used.
#$ -pe smp 2         # number of parallel processes to run (defines how many GPUs)
# -fin Example_Usage.py   # copy input file in the running directory
# -fin GeometryCheck.py
# -fin DataAnalysis.py 

module load anaconda3 # loads the appropriate terachem variables
export OMP_NUM_THREADS=2 # number of threads to be spawned (needs to be the same as -smp)

echo "running cage jobs!"
cd /home/manets12/DavidMimichrome 
# run the actual calculation and pipe output to empty_nanocage.out
conda run -n molsimp python Example_Usage.py > Efield_multiwfn_ExampleTest.log
