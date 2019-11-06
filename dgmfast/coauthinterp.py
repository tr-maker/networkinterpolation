"""
Interpolate with the snapshots from the coauthor datafiles. Save edge edits in "allgraphedits.npy".

Dependencies: snap[0-8].npz
"""
import numpy as np
import scipy.sparse as sp
import dgmfast as dgm
import time


numSnaps = 9

# Interpolate.
dtarget = 0
dtrigger = 0
slowness = 1
allgraphedits = []
snaps = [None] * numSnaps
for i in range(numSnaps):
    snaps[i] = sp.load_npz("snap" + str(i) + ".npz")

t0 = time.time()
for i in range(numSnaps-1):
    graphedits = dgm.dgmfast(snaps[i+1], snaps[i], dtarget, slowness, dtrigger)
    allgraphedits.append(graphedits)
    # print(i)
t1 = time.time()
print(t1 - t0)

np.save("allgraphedits", allgraphedits)

# ALTERNATIVELY:
# Interpolate.
# dtarget = 0
# dtrigger = 0
# slowness = 1
# allgraphedits = []
# currsnap = sp.load_npz("snap0.npz")
# for i in range(numSnaps-1):
#     nextsnap = sp.load_npz("snap" + str(i+1) + ".npz")
#     graphedits = dgm.dgmfast(nextsnap, currsnap, dtarget, slowness, dtrigger)
#     allgraphedits.append(graphedits)
#     currsnap = nextsnap
#     # print(i)
# np.save("allgraphedits", allgraphedits)
