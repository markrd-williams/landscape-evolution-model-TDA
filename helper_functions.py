#!/usr/bin/env python3

import os
import scipy.io
from scipy.spatial import distance_matrix
from scipy import sparse
from ripser import ripser
import numpy as np

def lower_star_distance_matrix(img):
    """
    Construct a lower star filtration on an image
    Parameters
    ----------
    img: ndarray (M, N)
        An array of single channel image data
    Returns
    -------
    I: ndarray (K, 2)
        A 0-dimensional persistence diagram corresponding to the sublevelset filtration
    """
    m, n = img.shape

    idxs = np.arange(m * n).reshape((m, n))

    I = idxs.flatten()
    J = idxs.flatten()
    V = img.flatten()

    # Connect 8 spatial neighbors
    tidxs = np.ones((m + 2, n + 2), dtype=np.int64) * np.nan
    tidxs[1:-1, 1:-1] = idxs

    tD = np.ones_like(tidxs) * np.nan
    tD[1:-1, 1:-1] = img

    for di in [-1, 0, 1]:
        for dj in [-1, 0, 1]:

            if di == 0 and dj == 0:
                continue

            thisJ = np.roll(np.roll(tidxs, di, axis=0), dj, axis=1)
            thisD = np.roll(np.roll(tD, di, axis=0), dj, axis=1)
            thisD = np.maximum(thisD, tD)

            # Deal with boundaries
            boundary = ~np.isnan(thisD)
            thisI = tidxs[boundary]
            thisJ = thisJ[boundary]
            thisD = thisD[boundary]


            I = np.concatenate((I, thisI.flatten()))
            J = np.concatenate((J, thisJ.flatten()))
            V = np.concatenate((V, thisD.flatten()))

    return sparse.coo_matrix((V, (I, J)), shape=(idxs.size, idxs.size))

def load_tops():
    directory = os.fsencode("./LEM/")
    data = {}
    for file in os.listdir(directory):
        if not file.endswith(b".nc"):
            continue
        filename = os.fsdecode(file[:-3])
        with scipy.io.netcdf_file(directory + file) as lem, scipy.io.netcdf_file(directory + b"/info/" + file[:-3] + b"params.nc") as params:
                top_list = lem.variables['topography__elevation'][:].copy()
                cell_size = float(params.variables['grid__length'][0])
                shape = top_list[0].shape

                data[filename] = {
                        "tops": top_list,
                        "cell_size": cell_size,
                        "dgms": [],
                        }

                # Augment dgms
                for top in top_list:
                    dist = lower_star_distance_matrix(-top)
                    dgm = ripser(dist, distance_matrix= True, do_cocycles=True) 
                    component_centres = np.array(dgm['cocycles'][0])[:, :, 0]
                    (x, y) = np.unravel_index(component_centres, shape)
                    x = x #* cell_size ignored for now as giving bad results
                    y = y #* cell_size
                    h0_barcode = dgm['dgms'][0]
                    data[filename]["dgms"].append( np.hstack((h0_barcode,) +(x,) + (y,)) )
    return data
