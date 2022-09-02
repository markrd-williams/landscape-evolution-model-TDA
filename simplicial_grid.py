import numpy as np
import itertools
from vineyard.vineyard import Simplex
from vineyard.vineyard import Vineyard

def grid(data):
    """
    data is out*n*m ndarray.
    represents n*m time series
    """

    (out,n,m) = data.shape

    # return the flattend index of t
    def i(t):
        return t[0] * n + t[1]

    # return the series of values of the vertex t
    def v(t):
        return data[:, t[0], t[1]]

    def simplex(coords):
        max_arr = np.full(out, -np.inf)
        for (x,y) in coords:
            data_arr = v((x,y))
            max_arr = np.maximum(max_arr, data_arr)
        max_arr = list(max_arr)
        return Simplex(list(map(i, coords)), max_arr)

    simplices = []
    for (x,y) in itertools.product(range(n), range(m)):
        simplices.append(simplex([(x,y)]))

        # these are coords of the potential simplices we need to add
        A = x + 1 < n and y - 1 >= 0
        B = x + 1 < n
        C = x + 1 < n and y + 1 < m
        D = y + 1 < m

        # Only add one dimensional simplices as only care about 0 homology
        if A:
            simplices.append(simplex([(x,y), (x+1, y-1)]))
        if B:
            simplices.append(simplex([(x,y), (x+1, y)]))
        if C:
            simplices.append(simplex([(x,y), (x+1, y+1)]))
        if D:
            simplices.append(simplex([(x,y), (x, y+1)]))
    return simplices

if __name__=='__main__':
    data = np.array( [[[-5, -9, 1], [1, 7, 2], [0,1, 3]], [[-5, -9, 1], [1, 7, 2], [0,1, 3]]] )
    simplices = grid(data)
    print(len(simplices))
    for s in simplices:
        print(s)

