import os
import sys
import subprocess as sp
from collections import defaultdict

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pyfastx import Fastq
#from kneed import KneeLocator

from .utils import read_config, readme_parser, dir_check, file_check



class PreFlight:
    def __init__(self, config):
        self.config = config
        self.sample = self.config['sample']
        proj_dir = self.config['outdir']

        self.data_dir = os.path.join(proj_dir, 'clean_data')
        dir_check(self.data_dir)

        self.res_dir = os.path.join(proj_dir, 'Result')
        dir_check(self.res_dir)

    def run(self):
        prededup = self.pre_dedup()
        fq_dedup = self.dedup(prededup)
        fq_valid = self.ready_fq(fq_dedup)
        return fq_valid

    def output_fq(self, cell_num):
        fq_dedup = os.path.join(self.data_dir, f'{self.sample}.dedup.fq')
        reads_ctx = os.path.join(self.data_dir, f'{self.sample}_read_counts.tsv')
        fumi_ctx = os.path.join(self.data_dir, f'{self.sample}_UMI_counts.tsv')
        prefix = os.path.join(self.data_dir, f'{self.sample}_{cell_num}')

        dir_name = os.path.abspath(os.path.dirname(__file__))
        reads_retriev = os.path.join(dir_name, 'readsRetriev2')

        sp.run([
            reads_retriev,
            '-bc', fumi_ctx,
            '-cells', str(cell_num),
            '-fq', fq_dedup,
            '-o', prefix,
            '-tab', 'off'
            ])

        self.violin_plot(fumi_ctx, cell_num, 'UMI')
        self.violin_plot(reads_ctx, cell_num, 'Reads')
        return None

    def clean_data(self):
        fmd = self.config['sampleInfo']
        fq1, fq2 = readme_parser(fmd)
        fqout = os.path.join(self.data_dir, f'{self.sample}_')

        anchoradp = self.config['anchoradp']

        cmd = f'{anchoradp} -1<(zcat {fq1}) -2<(zcat {fq2}) -o {fqout}'

        param = self.config['anchoradpPara']
        if param:
            cmd += f' {param}'

        return fqout

    def pre_dedup(self):
        read_counts = defaultdict(int)
        fq1 = os.path.join(self.data_dir, f'{self.sample}_1.fq')
        fq2 = os.path.join(self.data_dir, f'{self.sample}_2.fq')

        file_check(fq1)
        file_check(fq2)

        prededup = os.path.join(self.data_dir, f'{self.sample}.prededup.fq')
        fh = open(prededup, 'w')
        for item1, item2 in zip(Fastq(fq1, build_index=False, full_name=True), \
            Fastq(fq2, build_index=False, full_name=True)):
            id_1, seq_1, qual_1 = item1
            id_2, seq_2, qual_2 = item2

            #if len(seq_2) < 60:
            #    continue
            cb = seq_1[0:20]
            read_counts[cb] += 1

            new_seq = seq_1 + seq_2
            new_qual = qual_1 + qual_2

            fh.write(f'@{id_2}\n{new_seq}\n+\n{new_qual}\n')
        fh.close()

        reads_ctx = os.path.join(self.data_dir, f'{self.sample}_read_counts.tsv')
        
        x = []
        y = []
        n = 1
        with open(reads_ctx, 'w') as fh:
            for cb, val in sorted(read_counts.items(), key=lambda item: item[1], reverse=True):
                fh.write(f'{cb}\t{val}\n')
                y.append(val)
                x.append(n)
                n += 1
        self.scatter_plot(reads_ctx, 'reads')
        return prededup
        
    def dedup(self, prededup):
        fq_dedup = os.path.join(self.data_dir, f'{self.sample}.dedup.fq')

        nubeam_dedeup = self.config['nubeam_dedup']

        sp.run([
            nubeam_dedeup,
            '-i', prededup,
            '-o', fq_dedup
        ])

        return fq_dedup

    def ready_fq(self, fq_dedup, cell_num=50000):
        fq_final = os.path.join(self.data_dir, f'{self.sample}_final_2.fq')
        reads_ctx = os.path.join(self.data_dir, f'{self.sample}_read_counts.tsv')
        fumi_ctx = os.path.join(self.data_dir, f'{self.sample}_UMI_counts.tsv')
        prefix = os.path.join(self.data_dir, self.sample)

        dir_name = os.path.abspath(os.path.dirname(__file__))
        reads_retriev = os.path.join(dir_name, 'readsRetriev2_tmp')

        sp.run([
            reads_retriev,
            '-bc', reads_ctx,
            '-cells', str(cell_num),
            '-fq', fq_dedup,
            '-o', prefix,
            '-tab', 'on'
            ])

        self.scatter_plot(fumi_ctx, 'UMI')

        return fq_final

    def fq_validbc(self, cell_num):
        fumi_ctx = os.path.join(self.data_dir, f'{self.sample}_UMI_counts.tsv')
        fq_valid = os.path.join(self.data_dir, f'{self.sample}_top{cell_num}.fq')
        fq_final = os.path.join(self.data_dir, f'{self.sample}.final.fq')

        if not cell_num:
            cell_num = 1000

        cells = []
        n = 0
        with open(fumi_ctx) as fh:
            for line in fh:
                if n > cell_num:
                    break
                arr = line.split()
                cells.append(arr[0])
                n += 1
        
        with open(fq_valid, 'w') as fh:
            for item in Fastq(fq_final, build_index=False, full_name=True):
                id_, seq, qual = item
                _, cb = id_.split('_')
                if cb not in cells:
                    continue
                fh.write(f'@{id_}\n{seq}\n+\n{qual}\n')
        return fq_valid

    def scatter_plot(self, fdata, ylabel):
        df = pd.read_csv(fdata, header=None, sep='\t')
        df['Index'] = range(1, df.shape[0]+1)
        #kn = KneeLocator(x, y, curve='convex', direction='decreasing')
        #print(kn.knee)

        x = df.loc[:, 'Index']
        y = df.loc[:, 1]

        fig = plt.figure(figsize =(7, 8))
        ax = fig.subplots()
        ax.scatter(x, y, color='blue', s = 0.5)
        ax.set_title(self.sample)
        ax.set_ylabel(ylabel)
        ax.set_xlabel('Barcode')
        ax.loglog(base = 10)
        ylim = max(y)
        #if kn.knee:
        #    ax.vlines(kn.knee, 0, ylim, colors = 'blue', linestyles = 'dashed')
        fig.tight_layout()
        plt.savefig(os.path.join(self.data_dir, f'{self.sample}_{ylabel}_ScatterPlot.png'))

        return None

    def violin_plot(self, fdata, cell_num, ylabel):
        df = pd.read_csv(fdata, header=None, sep='\t')
        df = df.iloc[0:cell_num, :]
        values = df.loc[:, 1]

        fig = plt.figure(figsize =(7, 8))
        ax = fig.subplots()
        sns.violinplot( y= values, inner='box', saturation=10)
        ax.set_ylabel(ylabel)
        plt.savefig(os.path.join(self.res_dir, f'{self.sample}_{cell_num}_vioPlot.png'))



