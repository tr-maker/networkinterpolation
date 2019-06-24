% Extrapolate, using preferential attachment/detachment, with the snapshots 
% from the coauthor datafiles. Save mean and global clustering 
% coefficients, along with step indices, in 'coauthextrappref.mat'.
%
% Dependencies: coauthsnaps.mat, coauthinterp.mat, extrappref.m

load('coauthsnaps.mat','snaps')  % cell array of adjacency matrix snapshots
load('coauthinterp.mat','stepsi','period')
numSnaps = length(snaps);
numNodes = size(snaps{1},1);

% EXTRAPOLATE
% Preferential attachment/detachment
acp = cell(numSnaps,1);  % mean clustering coefficients of extrapolating graphs
gcp = cell(numSnaps,1);  % global clustering coefficients of extrapolating graphs
stepsp = cell(numSnaps,1);  % step indices
for i = 1:numSnaps-1
    [acp{i},gcp{i},stepsp{i}] = extrappref(snaps{i+1},snaps{i},period);
    % Kludge to not connect the graph between snapshots.
    acp{i}(end+1) = NaN;
    gcp{i}(end+1) = NaN;
    stepsp{i}(end+1) = stepsp{i}(end);
    disp(['Preferential: ' num2str(i) ' out of ' num2str(numSnaps-1) ' done'])
end
[acp{numSnaps},gcp{numSnaps},stepsp{numSnaps}] = extrappref(snaps{numSnaps},...
    snaps{numSnaps},period);
% Make the starts of the step indices align with the interpolation.
for i = 2:numSnaps
    stepsp{i} = stepsp{i} + stepsi{i}(1);
end
save('coauthextrappref.mat','acp','gcp','stepsp')