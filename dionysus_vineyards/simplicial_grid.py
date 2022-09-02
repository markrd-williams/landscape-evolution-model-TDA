import numpy as np
import itertools
import sys
import scipy.io

def grid(data):
    """
    data is out*n*m ndarray.
    represents n*m time series
    """

    (out,n,m) = data.shape

    # return the flattend index of t
    def i(x, y):
        return x * m + y

    simplices = []
    vertices = data.reshape(out, n*m)
    for (x,y) in itertools.product(range(n), range(m)):
        simplices.append([i(x,y)])

        # these are coords of the potential simplices we need to add
        A = x + 1 < n and y - 1 >= 0
        B = x + 1 < n
        C = x + 1 < n and y + 1 < m
        D = y + 1 < m

        # Only add one dimensional simplices as only care about 0 homology
        if A:
            simplices.append([i(x,y), i(x+1, y-1)])
        if B:
            simplices.append([i(x,y), i(x+1, y)])
        if C:
            simplices.append([i(x,y), i(x+1, y+1)])
        if D:
            simplices.append([i(x,y), i(x, y+1)])
    simplices.sort(key=len)
    return simplices, vertices 

if __name__=='__main__':
    filename = sys.argv[1]

    with scipy.io.netcdf_file(filename) as f:
        data = f.variables['topography__elevation'][100:, 60:140, 60:140].copy()

    simplices, values = grid(data)

    with open("complex", "w") as f:
        for s in simplices:
            f.write(" ".join(map(str, s)))
            f.write("\n")

    with open("values", "w") as f:
        for v in values:
            f.write(" ".join(map(str, v)))
            f.write("\n")
