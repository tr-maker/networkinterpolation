% Read in time-stamped edges from the datafile 'emaileucore.txt'.
% Create two snapshots: the empty initial starting graph and the 
% final accumulated graph.
% Save adjacency matrices of the snapshots in a cell array in 
% 'emaileucoresnaps.mat'.
%
% Dependencies: emaileucore.txt

% Read in edges from datafile.
edges = dlmread('emaileucore.txt');
numNodes = max(edges,[],'all')+1;

numSnaps = 2;
snaps = cell(1,numSnaps);
snaps{1} = sparse(numNodes,numNodes);
graph = sparse(edges(:,1)+1,edges(:,2)+1,1,numNodes,numNodes);
graph = graph .* ~eye(size(graph));  % zero out the diagonal entries
graph = graph + graph';
graph = spones(graph);
snaps{2} = graph;

save('emaileucoresnaps.mat','snaps')