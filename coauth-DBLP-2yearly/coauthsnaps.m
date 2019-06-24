% Read in directed edges from coauthor datafiles.
% Save adjacency matrices of snapshots in 'coauthsnaps.mat'.
%
% Dependencies: coauth[1-9].txt

numSnaps = 9;
snaps = cell(1,numSnaps);
% Read in snapshots from datafiles.
% Create snapshot adjacency matrices.
% Note: the adjacency matrices are undirected.
numNodes = 0;
for i = 1:numSnaps
    edges = dlmread(['coauth' num2str(i) '.txt']);
    if max(edges,[],'all') > numNodes
        numNodes = max(edges,[],'all');
    end
end
for i = 1:numSnaps
    edges = dlmread(['coauth' num2str(i) '.txt']);
    graph = sparse(edges(:,1),edges(:,2),1,numNodes,numNodes);
    graph = graph + graph';
    %disp(diag(graph))
    graph = spones(graph);
    snaps{i} = graph;
end
save('coauthsnaps.mat','snaps')