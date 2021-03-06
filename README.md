# networkinterpolation
Code and formatted datasets for the paper "Network Interpolation" ([SIAM online](https://epubs.siam.org/doi/abs/10.1137/19M1268380), [arXiv](https://arxiv.org/abs/1905.01253)) by Thomas Reeves, Anil Damle, and Austin R. Benson.

``dgmfast.py`` in the folder "dgmfast" runs the network interpolation model until a specified graph edit distance to the target graph is reached. It is the fastest implementation.

``dgmtrigger.m`` also runs the network interpolation model until a specified graph edit distance to the target graph is reached, but is not as fast for large graphs.

``dgmrun.m`` runs the network interpolation model for a fixed number of steps.

The program ``twotothreeblock.m`` contains the synthetic experiments described in Section 4 of the paper.

The folders "collegemsg", "emaileucore", "emaileucoretemporal", "vdbunt", and "coauth-DBLP-2yearly" contain the real-world experiments described in Section 5 of the paper.
