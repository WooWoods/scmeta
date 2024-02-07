import os
import sys
import argparse

from .utils import read_config
from .preflight import PreFlight
from .classify import Classifier, run_kraken
from .report import Report



def main():
    AP = argparse.ArgumentParser(
            description="microbiome scRNA seq analysis pipeline.",
            formatter_class=argparse.RawTextHelpFormatter,
            )
    AP.add_argument('--cfg', metavar='', help='config file, refer to example for details', required=True)
    AP.add_argument('--report', action='store_true', help='Re-analysis with a specific cell number')
    AP.add_argument('--cellnum', type=int, help='Number of cells')

    args = AP.parse_args()

    fcfg = args.cfg
    config = read_config(fcfg)

    obj_pre = PreFlight(config)

    if args.report:
        cell_num = args.cellnum
        if not cell_num:
            print("Error: cellnum not provided.\n")
            AP.print_help()
            sys.exit(1)
        obj_pre.output_fq(cell_num)
        reporter = Report(config)
        reporter.report()
        reporter.report2(cell_num)

    else:
        fq_valid = obj_pre.run()

        kreport, koutput = run_kraken(fq_valid, config)
        
        classifier = Classifier(config, kreport)
        classifier.worker_pool(koutput)

        reporter = Report(config)
        reporter.report()

