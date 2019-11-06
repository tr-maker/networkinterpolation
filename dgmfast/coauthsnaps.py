"""
Read in directed edges from coauthor datafiles.
Save strictly upper triangular adjacency matrices of snapshots in "snap[0-8].npz".

Dependencies: coauth[1-9].txt
"""
import scipy.sparse as sp

numSnaps = 9
# Read in snapshots from datafiles.
# Create snapshot adjacency matrices.
# Note: the adjacency matrices are undirected.
numNodes = 0
for i in range(numSnaps):
    f = open("coauth" + str(i + 1) + ".txt", "r")
    for line in f:
        edges = [int(t) for t in line.split(" ")]
        if max(edges) > numNodes:
            numNodes = max(edges)
for i in range(numSnaps):
    snap = sp.dok_matrix((numNodes, numNodes), dtype=int)  # Fast insertion
    f = open("coauth" + str(i + 1) + ".txt", "r")
    for line in f:
        edges = [int(t) for t in line.split(" ")]
        edges.sort()
        assert(edges[0] < edges[1])
        snap[edges[0] - 1, edges[1] - 1] = 1  # Minimum vertex in the dataset is 1, not 0
    snap = snap.tocsc()
    # print(snap)
    sp.save_npz("snap" + str(i), snap)
