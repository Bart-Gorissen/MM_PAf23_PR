#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:01:00

repeats=10 # number of repetitions
L_choices=(10) # inlink number
p_choices=(0.85) # follow probability
u_choices=(0 1) # random and e_1
D_choices=(0 1 2 3) # all options for communication (fixed u=0, pGx=2)
pGx_choices=(0 2) # all choices for pGx prodcut (fixed u=0, D=2)

ns=(10 100 1000)
ps=(1 2 4 8)

mkdir -p /scratch-shared/$USER
cd /scratch-shared/#USER
mkdir out

for n in "${ns[@]}"
do
    for p in "${ps[@]}"
    do
        echo $n > input.in
        echo ${L_choices[0]} >> input.in
        echo ${p_choices[0]} >> input.in
        echo $p >> input.in
        echo -8 >> input.in

        for i in `seq 1 $repeats`; do
            $HOME/projects/MM_PAf23_PR/bin/main -f=input.in
        done
    done
done

cp -r /out $SLURM_SUBMIT_DIR/output