import os
import sys
from collections import Counter, defaultdict

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from .matrix import CountMatrix
from .utils import read_config, readme_parser, dir_check


class Report:
    def __init__(self, config):
        self.config = config
        filter_thresh = config.get('filter_threshold', None)
        try:
            self.filter_tresh = int(filter_tresh.strip())
            self.do_filter = True
        except Exception:
            self.do_filter =False

        self.sample = self.config['sample']
        proj_dir = self.config['outdir']
        self.res_dir = os.path.join(proj_dir, 'Result')
        self.data_dir = os.path.join(proj_dir, 'clean_data')
        self.braken_dir = os.path.join(self.res_dir, 'braken_report')
        self.braken_dir_g = os.path.join(self.res_dir, 'braken_report_g')
        self.mat_dir = os.path.join(self.res_dir, 'matrix')
        dir_check(self.mat_dir)

    def report(self):
        fresult = os.path.join(self.res_dir, f'{self.sample}_sc_allot.result')
        freport = os.path.join(self.res_dir, f'{self.sample}_sc_taxonomy.report')
        freportg = os.path.join(self.res_dir, f'{self.sample}_sc_taxonomy.G.report')
        fh1 = open(fresult, 'w')
        fh2 = open(freport, 'w')
        fh3 = open(freportg, 'w')
        fh1.write("barcode\tname\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n")
        fh2.write("barcode\tname\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n")
        fh3.write("barcode\tname\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n")

        bc_counts = defaultdict(lambda : defaultdict(int))
        features = set()

        for bout in os.listdir(self.braken_dir):
            fbout = os.path.join(self.braken_dir, bout)
            if not os.path.exists(fbout):
                continue
            cb, _ = bout.split('.')
            items = []
            with open(fbout) as fh4:
                fh4.readline()
                for line in fh4:
                    arr = line.strip().split('\t')
                    feature = arr[0].replace(' ', '_')
                    bc_counts[cb][feature] += int(arr[5])
                    features.add(feature)

                    if self.do_filter:
                        if int(arr[5]) < self.filter_thresh:
                            continue
                        items.append(arr)
                    else:
                        items.append(arr)
            if not items:
                continue

            items = sorted(items, key=lambda x: x[-1], reverse=True)

            fh2.write(f'{cb}\t' + '\t'.join(items[0]) + '\n')

            for i, item in enumerate(items[0:5]):
                if i == 0:
                    fh1.write(f'{cb}\t' + '\t'.join(item) + '\n')
                else:
                    fh1.write('-\t' + '\t'.join(item) + '\n')
        fh1.close()
        fh2.close()

        outdir_raw = os.path.join(self.mat_dir, 'raw')
        dir_check(outdir_raw)
        cMatrix = CountMatrix(bc_counts, features)
        cMatrix.save_mex(outdir_raw)

        for bout in os.listdir(self.braken_dir_g):
            fbout = os.path.join(self.braken_dir_g, bout)
            if not os.path.exists(fbout):
                continue
            cb, *tmp = bout.split('.')
            items = []
            with open(fbout) as fh5:
                fh5.readline()
                for line in fh5:
                    arr = line.strip().split('\t')
                    if int(arr[5]) < 10:
                        continue
                    #if float(arr[-1]) < 0.05:
                    #    continue
                    items.append(arr)
            items = sorted(items, key=lambda x: x[-1], reverse=True)
            if not items:
                continue

            fh3.write(f'{cb}\t' + '\t'.join(items[0]) + '\n')
        fh3.close()

        self.report2()
        return freport

    def report2(self, cell_num='all'):
        freport = os.path.join(self.res_dir, f'{self.sample}_sc_taxonomy.report')
        freportg = os.path.join(self.res_dir, f'{self.sample}_sc_taxonomy.G.report')

        fumi_ctx = os.path.join(self.data_dir, f'{self.sample}_UMI_counts.tsv')
        reads_ctx = os.path.join(self.data_dir, f'{self.sample}_read_counts.tsv')

        df = pd.read_csv(freport, header=0, index_col=0, sep='\t')
        dfg = pd.read_csv(freportg, header=0, index_col=0, sep='\t')

        if cell_num == 'all':
            self.pie_plot(df, cell_num)
            self.violin_plot(df, cell_num)
            return

        cells = set()
        n = 0
        with open(fumi_ctx) as fh:
            for line in fh:
                if n > cell_num:
                    break
                arr = line.split()
                cells.add(arr[0])
                n += 1
        co_cells = cells & set(df.index)
        co_cellsg = cells & set(dfg.index)

        freport2 = os.path.join(self.res_dir, f'{self.sample}_sc_taxonomy_{cell_num}_bc.report')
        freportg2 = os.path.join(self.res_dir, f'{self.sample}_sc_taxonomy_G_{cell_num}_bc.report')

        df_cells = df.loc[list(co_cells), :]
        dfg_cells = dfg.loc[list(co_cellsg), :]

        df_cells.to_csv(freport2, sep='\t')
        dfg_cells.to_csv(freportg2, sep='\t')

        self.pie_plot(df_cells, cell_num)
        self.violin_plot(df_cells, cell_num)

    def pie_plot(self, df, cell_num):
        counts = Counter(df['name'].values)
        data = counts.values()
        names = counts.keys()

        fig = plt.figure(figsize =(13, 8))
        ax = fig.subplots()
        wedges, texts = ax.pie(data, textprops=dict(color="w"))
        ax.legend(wedges, names,
                title="Species",
                loc="center left",
                bbox_to_anchor=(1, 0, 0.5, 1))
        plt.tight_layout()
        plt.savefig(os.path.join(self.res_dir, f'{self.sample}_{cell_num}_piePlot.png'))


    def violin_plot(self, df, cell_num):
        fig = plt.figure(figsize =(7, 8))
        ax = fig.subplots()
        sns.violinplot( y= df['fraction_total_reads'], inner='box', saturation=10)
        ax.set_ylabel('Purity')
        plt.savefig(os.path.join(self.res_dir, f'{self.sample}_{cell_num}_vioPlot.png'))



