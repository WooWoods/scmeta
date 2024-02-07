import os
import sys
import gzip

import numpy as np
from scipy import sparse, io
from scipy.sparse import csr_matrix


class CountMatrix:
    def __init__(self, bc_counts, features):
        self.bc_counts = bc_counts
        self.barcodes = list(bc_counts.keys())
        self.features = features

        self.to_matrix()

    def to_matrix(self):
        self.rows = len(self.features)
        self.cols = len(self.barcodes)

        mat = np.zeros([self.rows, self.cols])

        for i, feature in enumerate(self.features):
            for j, bc in enumerate(self.barcodes):
                try:
                    counts = self.bc_counts[bc][feature]
                    mat[i, j] = counts
                except KeyError:
                    continue
        csr_mat = csr_matrix(mat)
        self.m = csr_mat.tocoo()

    def save_mex(self, outdir):

        metadata = '%%MatrixMarket matrix coordinate real general\n'
        out_matrix_fn = os.path.join(outdir, 'matrix.mtx.gz')
        out_barcodes_fn = os.path.join(outdir, 'barcodes.tsv.gz')
        out_features_fn = os.path.join(outdir, 'features.tsv.gz')

        with gzip.open(out_barcodes_fn, 'wb') as fh:
            for bc in self.barcodes:
                fh.write(f'{bc}\n'.encode())

        with gzip.open(out_features_fn, 'wb') as fh:
            for fet in self.features:
                fh.write(f'{fet}\n'.encode())

        with gzip.open(out_matrix_fn, 'wb') as fh:
            fh.write(np.compat.asbytes(metadata))
            fh.write(np.compat.asbytes("%\n"))
            fh.write(np.compat.asbytes('%i %i %i\n' % (self.rows, self.cols, self.m.nnz)))

            # write row, col, val in 1-based indexing
            for r, c, d in zip(self.m.row+1, self.m.col+1, self.m.data):
                fh.write(np.compat.asbytes(("%i %i %i\n" % (r, c, d))))



