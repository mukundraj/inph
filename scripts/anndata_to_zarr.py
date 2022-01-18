#!/home/mraj/anaconda3/envs/Rv4/bin/python

import sys
from anndata import read_h5ad
import zarr


arggs = str(sys.argv);
arggs = arggs.split(',')
print (arggs[2])
adata = read_h5ad(arggs[1].strip().strip("'"))

opfile = arggs[2].strip().strip("]").strip("'")
adata.write_zarr(opfile)


# Update the zarr file index

zfile = zarr.open(opfile, mode='r+') # no need to close explicitly - https://zarr.readthedocs.io/en/stable/tutorial.html
zfile.var._index[:] = zfile.var.gene_id[:]
