% Interpolate with the snapshots from the coauthor datafiles. Save 
% mean and global clustering coefficients, along with step indices and 
% sampling period, in 'coauthinterp.mat'.
%
% Dependencies: coauthsnaps.mat, dgm.m

load('coauthsnaps.mat','snaps')  % cell array of adjacency matrix snapshots
numSnaps = length(snaps);
numNodes = size(snaps{1},1);

% INTERPOLATE
dtarget = 0;  % target edit distance
slowness = 1;  % rate of approach to the target edit distance (higher values are slower)
dtrigger = 0;  % graph distance to the background graph that stops the dynamic graph model
period = 5000;  % sampling period (higher values require less storage)

aci = cell(numSnaps,1);  % mean clustering coefficient of interpolating graphs
gci = cell(numSnaps,1);  % global clustering coefficient of interpolating graphs
stepsi = cell(numSnaps,1);  % step indices
for i = 1:numSnaps-1
    [aci{i},gci{i},stepsi{i}] = dgm(snaps{i+1},snaps{i},dtarget,...
        slowness,dtrigger,period);
    disp(['Interpolation: ' num2str(i) ' out of ' num2str(numSnaps-1) ' done'])
end
[aci{numSnaps},gci{numSnaps},stepsi{numSnaps}] = dgm(snaps{numSnaps},...
    snaps{numSnaps},dtarget,slowness,dtrigger,period);
% Make the step indices cumulative.
for i = 2:numSnaps
    stepsi{i} = stepsi{i} + stepsi{i-1}(end);
end
save('coauthinterp.mat','aci','gci','stepsi','period')