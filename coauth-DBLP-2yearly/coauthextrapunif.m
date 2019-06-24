% Extrapolate, using uniform attachment/detachment, with the snapshots 
% from the coauthor datafiles. Save mean and global clustering 
% coefficients, along with step indices, in 'coauthextrapunif.mat'.
%
% Dependencies: coauthsnaps.mat, coauthinterp.mat, extrapunif.m

load('coauthsnaps.mat','snaps')  % cell array of adjacency matrix snapshots
load('coauthinterp.mat','stepsi','period')
numSnaps = length(snaps);
numNodes = size(snaps{1},1);

% EXTRAPOLATE
% Uniform attachment/detachment
acu = cell(numSnaps,1);  % mean clustering coefficients of extrapolating graphs
gcu = cell(numSnaps,1);  % global clustering coefficients of extrapolating graphs
stepsu = cell(numSnaps,1);  % step indices
for i = 1:numSnaps-1
    [acu{i},gcu{i},stepsu{i}] = extrapunif(snaps{i+1},snaps{i},period);
    % Kludge to not connect the graph between snapshots.
    acu{i}(end+1) = NaN;
    gcu{i}(end+1) = NaN;
    stepsu{i}(end+1) = stepsu{i}(end);
    disp(['Uniform: ' num2str(i) ' out of ' num2str(numSnaps-1) ' done'])
end
[acu{numSnaps},gcu{numSnaps},stepsu{numSnaps}] = extrapunif(snaps{numSnaps},...
    snaps{numSnaps},period);
% Make the starts of the step indices align with the interpolation.
for i = 2:numSnaps
    stepsu{i} = stepsu{i} + stepsi{i}(1);
end
save('coauthextrapunif.mat','acu','gcu','stepsu')