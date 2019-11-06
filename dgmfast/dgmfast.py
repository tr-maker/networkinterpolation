"""
Given an n by n undirected target graph and an n by n undirected starting graph, run our graph interpolation model until
the specified graph distance to the target graph is reached.
Return a list of graph edits made in order from the starting graph.
Each entry in the list is the integer i*n + j + 1 if vertex (i,j) is added to the graph
and the integer -i*n - j - 1 if vertex (i,j) is deleted from the graph.
b is the adjacency matrix of the target graph, assumed to be a CSC sparse matrix.
a is the adjacency matrix of the starting graph, assumed to be a CSC sparse matrix.
(Only the strictly upper triangular parts of the matrices are used.)
dtarget is the target edit distance between the target graph and the current graph.
slowness is the rate of approach to the target edit distance (higher values are slower).
dtrigger is the graph distance to the target graph that stops the dynamic graph model.

Notes: False edges are allowed.
This implementation is efficient, being (amortized) constant time per time step.
"""
import numpy as np
import scipy.sparse as sp
import randomDict as rd
import random as random


def dgmfast(b, a, dtarget, slowness, dtrigger):
    n = a.shape[0]

    # u is a RandomDict.
    # Each entry in u is i*n + j + 1 if (i,j) is not in the current graph but in the target graph
    # and -i*n - j - 1 if (i,j) is in the current graph but not in the target graph.
    umat = sp.triu(b-a, 1)
    # current number of advancing moves (edit distance between target graph and current graph)
    d = umat.nnz
    # current number of regressing moves
    duseless = n*(n-1)/2 - d
    u = rd.RandomDict()
    for i in range(d):
        iindex = umat.row[i]
        jindex = umat.col[i]
        if umat.data[i] == 1:
            u.add(iindex * n + jindex + 1)
        else:
            u.add(-iindex * n - jindex - 1)

    # graphedits is a list of graph edits in order.
    # Each entry in graphedits is i*n + j + 1 if (i,j) is added to the graph
    # and -i*n - j - 1 if (i,j) is deleted from the graph.
    graphedits = []

    if d <= dtrigger:
        return graphedits

    i = 1
    while True:
        # Make a random move according to the "Markov" chain.
        prb = np.nan
        if d == 0:
            prb = 0
        elif duseless == 0:
            prb = 1
        else:
            # Use a logistic function for d -> p.
            prb = 1/(1 + np.exp((-d + dtarget)/slowness))

        if np.random.uniform() < prb:
            # Do an advancing move.
            usefuledge = u.get_random()
            if usefuledge > 0:
                graphedits.append(usefuledge)
                u.remove(usefuledge)
                d = d - 1
                duseless = duseless + 1
            else:
                graphedits.append(usefuledge)
                u.remove(usefuledge)
                d = d - 1
                # The following happens because we allow false edges.
                duseless = duseless + 1
        else:
            # Do a regressing move (via rejection sampling).
            uselessedge = random.sample(range(n), 2)
            uselessedge.sort()
            while u.search(uselessedge[0] * n + uselessedge[1] + 1) is not None or \
                    u.search(-uselessedge[0] * n - uselessedge[1] - 1) is not None:
                uselessedge = random.sample(range(n), 2)
                uselessedge.sort()

            if b[uselessedge[0], uselessedge[1]] == 1:
                # Delete the edge.
                graphedits.append(-uselessedge[0] * n - uselessedge[1] - 1)
                u.add(uselessedge[0] * n + uselessedge[1] + 1)
                d = d + 1
                duseless = duseless - 1
            else:
                # Add the edge.
                # The following happens because we allow false edges.
                graphedits.append(uselessedge[0] * n + uselessedge[1] + 1)
                u.add(-uselessedge[0] * n - uselessedge[1] - 1)
                d = d + 1
                duseless = duseless - 1

        # print("Step")
        # print(i)
        # print(graphedits)
        # print()

        if d <= dtrigger:
            return graphedits

        i = i + 1


# Unit test, evolving from a clique on 3 nodes to a single edge 0--2 on 3 nodes.
row = np.array([0, 0, 1, 1, 2, 2])
col = np.array([1, 2, 0, 2, 0, 1])
data = np.array([1, 1, 1, 1, 1, 1])
a = sp.csc_matrix((data, (row, col)), shape=(3, 3))
row = np.array([0, 2])
col = np.array([2, 0])
data = np.array([1, 1])
b = sp.csc_matrix((data, (row, col)), shape=(3, 3))
dtarget = 0
dtrigger = 0
slowness = 1
dgmfast(b, a, dtarget, slowness, dtrigger)

"""
a = sp.dok_matrix((5, 5))
a[0, 0] = 1
a[0, 2] = 10
a[0, 2] = -10
a[0, 3] = 1234
a[0, 1] = 4
a.pop((0, 1))
print(a)
"""