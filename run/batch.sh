
for bin in ./bin/smic_100/gpuk20*; do echo $bin; done

# production
./batch.py ./bin/production_run_1000/smic_cpu 1 1 1 1 8






# components

export OMP_PROC_BIND=true
export OMP_NUM_THREADS=20
for bin in ./bin/smic_1000/cpu*; do ./batch.py $bin 20 20 1 1 5; done


export OMP_PROC_BIND=true
export OMP_NUM_THREADS=240
export MIC_OMP_PROC_BIND=true
export MIC_KMP_AFFINITY=SCATTER
for bin in ./bin/smic_100/mic*; do ./batch.py $bin 240 240 1 1 5; done




for bin in ./bin/smic_1000/gpu*; do ./batch.py $bin 14 14 1 1 5; done







# reps

export OMP_PROC_BIND=true
export OMP_NUM_THREADS=20
./batch.py ./bin/production_1000/smic_cpu* 1 80 1 1 8


export OMP_PROC_BIND=true
export OMP_NUM_THREADS=240
export MIC_OMP_PROC_BIND=true
export MIC_KMP_AFFINITY=SCATTER
./batch.py ./bin/production_1000/smic_mic* 1 960 1 1 8




./batch.py ./bin/production_1000/smic_gpu* 1 54 1 1 8

