% Read in time-stamped edges from the datafile 'collegemsg.txt'.
% Create two snapshots: the empty initial starting graph and the 
% accumulated graph 30 days later.
% Save adjacency matrices of the snapshots, and save the actual mean and 
% global clustering coefficients over time in 'collegemsgsnaps.mat'.
%
% Dependencies: collegemsg.txt, avgClusteringCoefficient.m

% Read in edges from datafile.
fileID = fopen('collegemsg.txt','r');
formatSpec = '%d %d %d';
sz = [3 Inf];
edges = fscanf(fileID,formatSpec,sz);
times = datetime(edges(3,:),'ConvertFrom','posixtime'); % timestamps
edges = edges(1:2,:); % edges in chronological order
numNodes = max(edges,[],'all');

% Create snapshot adjacency matrices.
% Note: the adjacency matrices are undirected.
timeWindow = duration(24*30,0,0);  % time window (duration(h,mi,s))
numTimes = 2;
while times(numTimes) - times(1) < timeWindow
    numTimes = numTimes + 1;
end
%%disp(numTimes)

numSnaps = 2;
snaps = cell(1,numSnaps);
snaps{1} = zeros(numNodes,numNodes);
snapfinal = zeros(numNodes,numNodes);

period = 50;
acreal = zeros(1,ceil(numTimes/period)+1);
ccreal = zeros(numNodes,ceil(numTimes/period)+1);
gcreal = zeros(1,ceil(numTimes/period)+1);
[acreal(1), ccreal(:,1), gcreal(1)] = avgClusteringCoefficient(snaps{1});
for i = 1:numTimes
    snapfinal(edges(1,i),edges(2,i)) = 1;
    snapfinal(edges(2,i),edges(1,i)) = 1;
    if mod(i,period) == 0
        disp(i)
        [acreal(i/period+1), ccreal(:,i/period+1), gcreal(i/period+1)] = avgClusteringCoefficient(snapfinal);
    end
end
if mod(numTimes,period) ~= 0
    [acreal(end), ccreal(:,end), gcreal(end)] = avgClusteringCoefficient(snapfinal);
end
snaps{2} = snapfinal;
save('collegemsgsnaps.mat','snaps','acreal','ccreal','gcreal','period','numTimes')