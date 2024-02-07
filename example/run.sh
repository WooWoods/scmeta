#!/bin/bash
#SBATCH -J scMeta
#SBATCH  -p cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10


source /public/home/wangycgroup/public/01_Pipeline/scMeta/venv/bin/activate
export LD_LIBRARY_PATH=/public/home/wangycgroup/public/software/lib2

if [ ! -d clean_data ]
then
	mkdir clean_data
fi

/public/home/wangycgroup/public/software/anchoradp.o -1<(zcat raw/YA-del_R1.fq.gz) -2<(zcat raw/YA-del_R2.fq.gz) -o clean_data/YA-del_ -p -a 60

python /public/home/wangycgroup/public/01_Pipeline/scMeta/scMeta.py --cfg config.ini
