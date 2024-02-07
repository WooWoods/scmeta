import os
import sys
import copy

import multiprocessing
import subprocess as sp
from collections import defaultdict
from multiprocessing import Pool, Process, Queue

from .utils import read_config, readme_parser, dir_check


def run_kraken(fq_valid, config):
    sample = config['sample']
    proj_dir = config['outdir']
    res_dir = os.path.join(proj_dir, 'Result')
    kraken2 = config['kraken']
    kdb = config['krakenDb']

    fq_classify = os.path.join(res_dir, 'output_classified.fq')
    fq_unclassify = os.path.join(res_dir, 'output_unclassified.fq')

    output = os.path.join(res_dir, f'{sample}_kraken.output')
    report = os.path.join(res_dir, f'{sample}_kraken.report')

    sp.run([
        kraken2,
        '--threads', config.get('process', '4'),
        '--db', kdb,
        '--unclassified-out', fq_unclassify,
        '--classified-out', fq_classify,
        '--report', report,
        '--output', output,
        fq_valid
    ])
    return report, output


class Classifier:
    def __init__(self, config, kreport):
        self.config = config
        self.p = int(self.config['process'])
        self.sample = self.config['sample']
        proj_dir = self.config['outdir']
        self.res_dir = os.path.join(proj_dir, 'Result')
        self.braken_dir = os.path.join(self.res_dir, 'braken_report')
        self.braken_dir_g = os.path.join(self.res_dir, 'braken_report_g')
        self.tmp_dir = os.path.join(self.res_dir, 'braken_tmp')
        dir_check(self.tmp_dir)
        dir_check(self.res_dir)
        dir_check(self.braken_dir)
        dir_check(self.braken_dir_g)
        self.braken = self.config['braken']
        self.kdb = self.config['krakenDb']
        self.keys = ['ratio', 'sum', 'counts', 'class', 'tax_id', 'name']

        self.taxan = Taxanomy(kreport)
    
    def cell_report(self, koutput):
        self.tmp_dir = tmp_dir
        dir_check(tmp_dir)

        data = defaultdict(lambda : defaultdict(int))
        with open(koutput) as fh:
            for line in fh:
                arr = line.split()
                _, cb = arr[1].split('_')
                data[cb][arr[2]] += 1

        for cb, counts in data.items():
            cb_report = os.path.join(tmp_dir, f'{cb}.kreport')
            result, total = self.cell_taxan(counts)
            with open(cb_report, 'w') as fh:
                for item in result:
                    arr = [str(item[k]) for k in keys]
                    fh.write('\t'.join(arr) + '\n')
            self.run_braken(cb_report, cb)

    def producer(self, bc_counts, bc_queue):
        for bc in bc_counts:
            bc_queue.put(bc)
        return None

    def consumer(self, bc_counts, bc_queue):
        #curr_proc = multiprocessing.current_process()
        while True:
            cb = bc_queue.get()
            if cb is None:
                break

            counts = bc_counts[cb]
            cb_report = os.path.join(self.tmp_dir, f'{cb}.kreport')
            result, total = self.cell_taxan(counts)
            with open(cb_report, 'w') as fh:
                for item in result:
                    arr = [str(item[k]) for k in self.keys]
                    fh.write('\t'.join(arr) + '\n')
            self.run_braken(cb_report, cb)

    def worker_pool(self, koutput):
        bc_queue = Queue(maxsize = 50000)

        bc_counts = defaultdict(lambda : defaultdict(int))
        with open(koutput) as fh:
            for line in fh:
                arr = line.split()
                _, cb = arr[1].split('_')
                bc_counts[cb][arr[2]] += 1

        p1 = Process(target=self.producer, args=(bc_counts, bc_queue, ))
        p1.start()

        jobs = []
        for i in range(self.p):
            p = Process(target=self.consumer, args=(bc_counts, bc_queue, ))
            p.start()
            jobs.append(p)

        p1.join()
        for i in range(self.p):
            bc_queue.put(None)
        for p in jobs:
            p.join()

        return None

    def tax_scan(self, taxon_tree, tax_id, counts):
        idx = 0
        tax_cts = 0
        for i, item in enumerate(taxon_tree):
            if item['tax_id'] == tax_id:
                idx = i
                item['hit'] = 1
                tax_cts = counts[tax_id]
                item['counts'] = tax_cts
                item['sum'] += tax_cts
                break

        def recur_search(taxon_tree, i, tax_cts):
            item = taxon_tree[i]

            if item['indent'] == 0:
                item['hit'] = 1
                return 

            j = item['father']
            if j:
                taxon_tree[j]['hit'] = 1
                taxon_tree[j]['sum'] += tax_cts
            else:
                indent = item['indent']

                j = i - 1
                while j >= 0:
                    item_prev = taxon_tree[j]
                    if item_prev['indent'] < indent:
                        taxon_tree[i]['father'] = j
                        taxon_tree[j]['hit'] = 1
                        taxon_tree[j]['sum'] += tax_cts
                        break
                    j -= 1
            return recur_search(taxon_tree, j, tax_cts)
        recur_search(taxon_tree, idx, tax_cts)

        return None

    def cell_taxan(self, counts):
        taxon_tree = copy.deepcopy(self.taxan.data)

        total = 0
        for tax_id in counts:
            total += counts[tax_id]
            self.tax_scan(taxon_tree, tax_id, counts)

        result = []
        for item in taxon_tree:
            if item['hit'] == 0:
                continue
            ratio = item['sum'] * 100/total
            item['ratio'] = f'{ratio:.2f}'
            item['name'] = ' ' * item['indent'] + item['name']
            result.append(item)
        
        result[0]['sum'] = result[0]['counts']
        ratio = result[0]['counts'] * 100/total
        result[0]['ratio'] = f'{ratio:.2f}'

        return result, total

    def run_braken(self, report, cb):
        output = os.path.join(self.braken_dir, f'{cb}.braken')
        output_g = os.path.join(self.braken_dir_g, f'{cb}.G.braken')
        b_log = os.path.join(self.tmp_dir, f'{cb}.braken.log')
        g_log = os.path.join(self.tmp_dir, f'{cb}.G.braken.log')
        fh = open(b_log, 'w')
        fhg = open(g_log, 'w')

        sp.run([
            self.braken, 
            '-d', self.kdb,
            '-i', report,
            '-o', output,
            '-r', '100'
        ], stdout=fh)

        sp.run([
            self.braken, 
            '-d', self.kdb,
            '-i', report,
            '-o', output_g,
            '-r', '100',
            '-l', 'G'
        ], stdout=fhg)

class Taxanomy:
    def __init__(self, kreport):
        self.data = []
        self.create_db(kreport)

    def create_db(self, kreport):
        with open(kreport) as fh:
            for line in fh:
                arr = line.strip().split('\t')
                item = self.create_item(arr)
                self.data.append(item)

    def search(self, tax_id):
        result = []
        def recur_search(data, i, result):
            item = data[i]

            if item['indent'] == 0:
                if result[-1]['name'] == 'root':
                    return result[0:-1]
                return result

            j = item['father']
            if j:
                father = data[j]
                result.append(father)
            else:
                indent = item['indent']

                j = i - 1
                while j > 1:
                    item_prev = data[j]
                    if item_prev['indent'] < indent:
                        data[i]['father'] = j
                        result.append(item_prev)
                        break
                    j -= 1
            return recur_search(data, j, result)

        idx = 0
        for i, item in enumerate(self.data):
            if item['tax_id'] == tax_id:
                idx = i
                break
        result.append(self.data[idx])
        recur_search(self.data, idx, result)
        return result

    def create_item(self, arr):
        sci_name = arr[-1].strip()
        tmp = arr[-1].split(' ')
        
        n = 0
        for i in tmp:
            if not i:
                n += 1

        item = {
                'tax_id': arr[4],
                'class' : arr[3],
                'name'  : sci_name,
                'indent': n,
                'father': '',
                'counts': 0, 
                'sum'   : 0, 
                'hit'   : 0
                }
        return item



if __name__ == '__main__':
    f = sys.argv[1]
    taxan = Taxanomy(f)
    result = taxan.search('83627')
    print(result)

