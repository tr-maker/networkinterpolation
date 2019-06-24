% Read in directed edges from van de Bunt social network datafiles.
% Save adjacency matrices of snapshots as a cell array in
% 'vdbuntsnaps.mat'.
%
% Dependencies: VRND32T[0-6].DAT

numSnaps = 7;
snaps = cell(1,numSnaps);
% Read in snapshots from datafiles.
% Create snapshot adjacency matrices.
% Note: the adjacency matrices are undirected.
for i = 1:numSnaps
    s = dlmread(['VRND32T' num2str(i-1) '.DAT']);
    s = s .* ~eye(size(s));  % zero out the diagonal entries
    
    % Create a friendship network.
    s(s == 2) = 1;
    s(s ~= 1) = 0;
    
    snaps{i} = double(s | s');
end
save('vdbuntsnaps.mat','snaps')