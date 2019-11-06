% Interpolate between the snapshots from the emaileucoretemporal datafile.
% Extrapolate to from the emaileucoretemproal datafile
% using the uniform attachment model, preferential attachment model, and 
% the triangle closing model.
% Plot mean and global clustering coefficients of real data, interpolation, 
% and extrapolations.
% Save the figure with the mean clustering coefficients as 
% 'emaileucoretemporalac.fig' and 'emaileucoretemporalac.eps'.
% Save the figure with the global clustering coefficients as 
% 'emaileucoretemporalgc.fig' and 'emaileucoretemporalgc.eps'.
% Save the experiment data in 'emaileucoretemporalanalysis.mat.'
% 
% Dependencies: emaileucoresnaps.mat, hittingtime.m, dgm.m, extrapunif.m,
% extrappref.m, extraptri.m

% cell array of adjacency matrix snapshots
load('emaileucoresnaps.mat','snaps','acreal','gcreal','stepsreal')
numSnaps = length(snaps);
numNodes = size(snaps{end},1);

% INTERPOLATE
aci = cell(numSnaps,1);
gci = cell(numSnaps,1);
stepsi = cell(numSnaps,1);
dtarget = 0;  % target edit distance
% For each interpolation, calculate best rate parameter (to the nearest
% multiple of 50) given the actual number of steps taken and the edit 
% distances betwen the snapshots.
slownesses = zeros(1,numSnaps-1);
for i = 1:numSnaps-1
    d0 = nnz(triu(snaps{i+1}-snaps{i},1)); % edit distance between snapshots
    if d0 == 0
        slownesses(i) = 0;
        continue;
    end
    hreal = stepsreal{i+1}(1) - stepsreal{i}(1); % number of steps taken between snapshots
    s = 0;
    while true
        [h1,~] = hittingtime(d0,0,s);
        if h1 >= hreal
            break;
        end
        s = s + 50;
    end
    if s == 0
        slownesses(i) = 0;
        continue;
    end
    [h2,~] = hittingtime(d0,0,s-50);
    if hreal - h2 <= h1 - hreal
        slownesses(i) = s-50;
    else
        slownesses(i) = s;
    end    
end
dtrigger = 0;  % graph distance to the target graph that stops the dynamic graph model
period = 100;  % sampling period (higher values require less storage)
for i = 1:numSnaps-1
    [aci{i},~,gci{i},stepsi{i}] = dgm(snaps{i+1},snaps{i},dtarget,...
        slownesses(i),dtrigger,period);
    % Kludge to not connect the graph between snapshots.
    aci{i}(end+1) = NaN;
    gci{i}(end+1) = NaN;
    stepsi{i}(end+1) = stepsi{i}(end);
end
[aci{numSnaps},~,gci{numSnaps},stepsi{numSnaps}] = ...
    dgm(snaps{numSnaps},snaps{numSnaps},dtarget,0,dtrigger,period);
% Make the starts of the step indices align with the real data.
for i = 2:numSnaps
    stepsi{i} = stepsi{i} + stepsreal{i}(1);
end
disp('Interpolation done')

% EXTRAPOLATE
% Uniform attachment/detachment
acu = cell(numSnaps,1);  % clustering coefficients of extrapolating graphs
gcu = cell(numSnaps,1);
stepsu = cell(numSnaps,1);  % step indices
% Start with a tiny clique for technical reasons.
startingclique = snaps{1};
startingclique(1:2,1:2) = ~eye(2);
[acu{1},gcu{1},stepsu{1}] = extrapunif(snaps{2},startingclique,period);
% Kludge to not connect the graph between snapshots.
acu{1}(end+1) = NaN;
gcu{1}(end+1) = NaN;
stepsu{1}(end+1) = stepsu{1}(end);
% Do the rest of the extrapolations.
for i = 2:numSnaps-1
    [acu{i},gcu{i},stepsu{i}] = extrapunif(snaps{i+1},snaps{i},period);
    % Kludge to not connect the graph between snapshots.
    acu{i}(end+1) = NaN;
    gcu{i}(end+1) = NaN;
    stepsu{i}(end+1) = stepsu{i}(end);
end
[acu{numSnaps},gcu{numSnaps},stepsu{numSnaps}] = extrapunif(snaps{numSnaps},...
    snaps{numSnaps},period);
% Make the starts of the step indices align with the real data.
for i = 2:numSnaps
    stepsu{i} = stepsu{i} + stepsreal{i}(1);
end
disp('Uniform done')

% Preferential attachment/detachment
acp = cell(numSnaps,1);  % clustering coefficients of extrapolating graphs
gcp = cell(numSnaps,1);
stepsp = cell(numSnaps,1);  % step indices
% Start with a tiny clique for technical reasons.
[acp{1},gcp{1},stepsp{1}] = extrappref(snaps{2},startingclique,period);
% Kludge to not connect the graph between snapshots.
acp{1}(end+1) = NaN;
gcp{1}(end+1) = NaN;
stepsp{1}(end+1) = stepsp{1}(end);
% Do the rest of the extrapolations.
for i = 2:numSnaps-1
    [acp{i},gcp{i},stepsp{i}] = extrappref(snaps{i+1},snaps{i},period);
    % Kludge to not connect the graph between snapshots.
    acp{i}(end+1) = NaN;
    gcp{i}(end+1) = NaN;
    stepsp{i}(end+1) = stepsp{i}(end);
end
[acp{numSnaps},gcp{numSnaps},stepsp{numSnaps}] = extrappref(snaps{numSnaps},...
    snaps{numSnaps},period);
% Make the starts of the step indices align with the real data.
for i = 2:numSnaps
    stepsp{i} = stepsp{i} + stepsreal{i}(1);
end
disp('Preferential done')

% Triangle closing
act = cell(numSnaps,1);  % clustering coefficients of extrapolating graphs
gct = cell(numSnaps,1);
stepst = cell(numSnaps,1);  % step indices
% Start with a tiny clique for technical reasons.
[act{1},gct{1},stepst{1}] = extraptri(snaps{2},startingclique,period);
% Kludge to not connect the graph between snapshots.
act{1}(end+1) = NaN;
gct{1}(end+1) = NaN;
stepst{1}(end+1) = stepsp{1}(end);
% Do the rest of the extrapolations.
for i = 2:numSnaps-1
    [act{i},gct{i},stepst{i}] = extraptri(snaps{i+1},snaps{i},period);
    % Kludge to not connect the graph between snapshots.
    act{i}(end+1) = NaN;
    gct{i}(end+1) = NaN;
    stepst{i}(end+1) = stepst{i}(end);
end
[act{numSnaps},gct{numSnaps},stepst{numSnaps}] = extraptri(snaps{numSnaps},...
    snaps{numSnaps},period);
% Make the starts of the step indices align with the real data.
for i = 2:numSnaps
    stepst{i} = stepst{i} + stepsreal{i}(1);
end
disp('Triangle closing done')

% SAVE
save('emaileucoretemporalanalysis.mat','stepsi','aci','gci',...
    'acu','gcu',...
    'acp','gcp',...
    'act','gct',...
    'period')

% PLOT
figure
hold on
p0 = plot(horzcat(stepsreal{:}),horzcat(acreal{:}),'Color',[0.4660    0.6740    0.1880],'LineWidth',2);
p1 = plot(horzcat(stepsi{:}),horzcat(aci{:}),'Color',[0    0.4470    0.7410],'LineWidth',1.25);
p2 = plot(horzcat(stepsu{:}),horzcat(acu{:}),'--','Color',[0.8500    0.3250    0.0980],'LineWidth',1);
p3 = plot(horzcat(stepsp{:}),horzcat(acp{:}),'--','Color',[0.9290    0.6940    0.1250],'LineWidth',1);
p4 = plot(horzcat(stepst{:}),horzcat(act{:}),'--','Color',[0.4940    0.1840    0.5560],'LineWidth',1);
for i = 1:numSnaps
    scatter(stepsreal{i}(1),acreal{i}(1),'r')
end
legend([p0 p1 p2 p3 p4],{'Real data','Interpolation','Extrapolation (uniform)',...
    'Extrapolation (preferential)','Extrapolation (triangle closing)'},...
    'Location','best')
title('Institution messages')
xlabel('Steps')
ylabel('Mean clustering coefficient')
set(gca,'fontsize',20)
saveas(gcf,'emaileucoretemporalac.eps','epsc')
saveas(gcf,'emaileucoretemporalac.fig')

figure
hold on
p0 = plot(horzcat(stepsreal{:}),horzcat(gcreal{:}),'Color',[0.4660    0.6740    0.1880],'LineWidth',2);
p1 = plot(horzcat(stepsi{:}),horzcat(gci{:}),'Color',[0    0.4470    0.7410],'LineWidth',1.25);
p2 = plot(horzcat(stepsu{:}),horzcat(gcu{:}),'--','Color',[0.8500    0.3250    0.0980],'LineWidth',1);
p3 = plot(horzcat(stepsp{:}),horzcat(gcp{:}),'--','Color',[0.9290    0.6940    0.1250],'LineWidth',1);
p4 = plot(horzcat(stepst{:}),horzcat(gct{:}),'--','Color',[0.4940    0.1840    0.5560],'LineWidth',1);
for i = 1:numSnaps
    scatter(stepsreal{i}(1),gcreal{i}(1),'r')
end
legend([p0 p1 p2 p3 p4],{'Real data','Interpolation','Extrapolation (uniform)',...
    'Extrapolation (preferential)','Extrapolation (triangle closing)'},...
    'Location','best')
title('Institution messages')
xlabel('Steps')
ylabel('Global clustering coefficient')
set(gca,'fontsize',20)
saveas(gcf,'emaileucoretemporalgc.eps','epsc')
saveas(gcf,'emaileucoretemporalgc.fig')