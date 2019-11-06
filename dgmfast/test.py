import numpy as np
import scipy.sparse as sp
import random as random

# a = sp.csc_matrix((5, 5))
# a[0, 0] = 1
# a[0, 2] = 10
# a[0, 2] = -10
# a[1, 1] = 4
# print(a)
# indices = a.nonzero()
# print(indices)
# indicesflat = np.ravel_multi_index(indices, (5, 5))
# print(indicesflat)
# print(np.random.uniform())
# r = random.sample(range(10), 2)
# r.sort()
# print(r)
row = np.array([0, 3, 1, 0])
col = np.array([0, 3, 1, 2])
# data = np.array([4, 5, 7, 9])
# m = sp.coo_matrix((data, (row, col)), shape=(4, 4))
# print(m.row)
# print(m.col)
# print(m.data)
# print(m.nnz)
# a = []
# a.append(3)
# s = "3 4"
# s = [int(t) for t in s.split(" ")]
# print(s[0])
# print(s[1])