% Read in time-stamped edges from the datafile 'emaileucoretemporal.txt'.
% Create 10 cumulative snapshots (roughly every 100 days in the dataset).
% Record mean and global clustering coefficients for the actual data
% every 100 timestamps.
% Save adjacency matrices of the snapshots and mean and global clustering 
% coefficients for the actual data in cell arrays in 
% 'emaileucoresnaps.mat'.
%
% Dependencies: emaileucoretemporal.txt, avgClusteringCoefficient.m

% Read in edges from datafile.
fileID = fopen('emaileucoretemporal.txt','r');
formatSpec = '%d %d %d';
sz = [3 Inf];
edges = fscanf(fileID,formatSpec,sz);
times = datetime(edges(3,:),'ConvertFrom','posixtime'); % timestamps
[times,perm] = sort(times);
edges = edges(1:2,perm) + 1; % edges in chronological order
numNodes = max(edges,[],'all');
numTimes = length(times);

numSnaps = 10;
snaps = cell(1,numSnaps);
snaps{1} = sparse(numNodes,numNodes);
% Create a growing network, with snapshots corresponding to
% roughly every 100 days.
index = 1;
graph = zeros(numNodes,numNodes);
timeIndices = zeros(1,numSnaps);
timeIndices(1) = 0;
for i = 1:numTimes
    graph(edges(1,i),edges(2,i)) = 1;
    graph(edges(2,i),edges(1,i)) = 1;
    if times(i) - times(1) > duration(index*100*24,0,0)
        % Save the snapshot.
        snap = sparse(graph);
        snap = snap + snap';
        snap = snap .* ~eye(size(snap));  % zero out the diagonal entries
        %disp(diag(snap))
        snap = spones(snap);
        snaps{index+1} = snap;
        timeIndices(index+1) = i;
        index = index + 1;
        disp(i)
    end
end
snap = sparse(graph);
snap = snap + snap';
snap = snap .* ~eye(size(snap));  % zero out the diagonal entries
%disp(diag(snap))
snap = spones(snap);
snaps{index+1} = snap;
timeIndices(index+1) = numTimes;

% Save statistics of the real graph.
% First record statistics for every 100 edits.
graph = zeros(numNodes,numNodes);
period = 100;
acrealarray = zeros(1,ceil(numTimes/period)+1); 
gcrealarray = zeros(1,ceil(numTimes/period)+1); 
stepsrealarray = zeros(1,ceil(numTimes/period)+1);
for i = 1:numTimes
    graph(edges(1,i),edges(2,i)) = 1;
    graph(edges(2,i),edges(1,i)) = 1;
    if mod(i,period) == 0
        snap = sparse(graph);
        snap = snap + snap';
        snap = snap .* ~eye(size(snap));  % zero out the diagonal entries
        %disp(diag(snap))
        snap = spones(snap);
        [acrealarray(i/period+1),~,gcrealarray(i/period+1)] = avgClusteringCoefficient(snap);
        stepsrealarray(i/period+1) = i;
        disp(i)
    end
end
if mod(numTimes,period) ~= 0
    snap = sparse(graph);
    snap = snap + snap';
    snap = snap .* ~eye(size(snap));  % zero out the diagonal entries
    %disp(diag(snap))
    snap = spones(snap);
    [acrealarray(end),~,gcrealarray(end)] = avgClusteringCoefficient(snap);
    stepsrealarray(end) = i;
end

% Then align the statistics records with the snapshots.
acreal = cell(numSnaps,1);
gcreal = cell(numSnaps,1);
stepsreal = cell(numSnaps,1);
index = 1;
for i = 1:ceil(numTimes/period)+1
    if stepsrealarray(i) >= timeIndices(index+1)
        index = index+1;
    end
    acreal{index} = [acreal{index} acrealarray(i)];
    gcreal{index} = [gcreal{index} gcrealarray(i)];
    stepsreal{index} = [stepsreal{index} stepsrealarray(i)];
end
for i = 1:numSnaps-1
    [acr,~,gcr] = avgClusteringCoefficient(snaps{i});
    acreal{i} = [acr acreal{i}];
    gcreal{i} = [gcr gcreal{i}];
    stepsreal{i} = [timeIndices(i) stepsreal{i}];
    [acr,~,gcr] = avgClusteringCoefficient(snaps{i+1});
    acreal{i} = [acreal{i} acr];
    gcreal{i} = [gcreal{i} gcr];
    stepsreal{i} = [stepsreal{i} timeIndices(i+1)];
end

save('emaileucoretemporalsnaps.mat','snaps','acreal','gcreal','stepsreal')