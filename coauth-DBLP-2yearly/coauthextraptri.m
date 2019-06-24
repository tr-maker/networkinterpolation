% Extrapolate, using triangle closing, with the snapshots from the coauthor 
% datafiles. Save mean and global clustering coefficients, along with step 
% indices, in 'coauthextraptri.mat'.
%
% Dependencies: coauthsnaps.mat, coauthinterp.mat, extraptri.m

load('coauthsnaps.mat','snaps')  % cell array of adjacency matrix snapshots
load('coauthinterp.mat','stepsi','period')
numSnaps = length(snaps);
numNodes = size(snaps{1},1);

% EXTRAPOLATE
% Triangle closing
act = cell(numSnaps,1);  % mean clustering coefficients of extrapolating graphs
gct = cell(numSnaps,1);  % global clustering coefficients of extrapolating graphs
stepst = cell(numSnaps,1);  % step indices
for i = 1:numSnaps-1
    [act{i},gct{i},stepst{i}] = extraptri(snaps{i+1},snaps{i},period);
    % Kludge to not connect the graph between snapshots.
    act{i}(end+1) = NaN;
    gct{i}(end+1) = NaN;
    stepst{i}(end+1) = stepst{i}(end);
    disp(['Triangle closing: ' num2str(i) ' out of ' num2str(numSnaps-1) ' done'])
end
[act{numSnaps},gct{numSnaps},stepst{numSnaps}] = extraptri(snaps{numSnaps},...
    snaps{numSnaps},period);
% Make the starts of the step indices align with the interpolation.
for i = 2:numSnaps
    stepst{i} = stepst{i} + stepsi{i}(1);
end
save('coauthextraptri.mat','act','gct','stepst')