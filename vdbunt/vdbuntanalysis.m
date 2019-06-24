% Interpolate and extrapolate (using the uniform attachment model, 
% preferential attachment model, and the triangle closing model) with the 
% snapshots from the van de Bunt social network datafiles.
% Draw four interpolating graphs equally spaced between the second-to-last 
% snapshot and the last snapshot (including the second-to-last snapshot and
% the last snapshot). Save these as 'vdbuntgraph[1-4].fig' and 
% 'vdbuntgraph[1-4].eos'.
% Plot mean and global clustering coefficients of all the .
% Save the figure with the mean clustering coefficients as 'vdbuntac.fig' 
% and 'vdbuntac.eps'.
% Save the figure with the global clustering coefficients as 'vdbuntgc.fig'
% and 'vdbuntgc.eps'.
% Save the experiment data in 'vdbuntanalysis.mat.'
%
% Dependencies: vdbuntsnaps.mat, dgmgraphs.m, avgClusteringCoefficient.m,
% extrapunif.m, extrappref.m, extraptricc.m

load('vdbuntsnaps.mat','snaps')  % cell array of adjacency matrix snapshots
numSnaps = length(snaps);
numNodes = size(snaps{1},1);

% INTERPOLATE
dtarget = 0;  % target edit distance
slowness = 1;  % rate of approach to the target edit distance (higher values are slower)
dtrigger = 0;  % graph distance to the background graph that stops the dynamic graph model
period = 1;  % sampling period (higher values require less storage).

% Get interpolating graphs.
graphsi = cell(numSnaps,1);  % interpolating graphs
stepsi = cell(numSnaps,1);  % step indices
for i = 1:numSnaps-1
    [graphsi{i},stepsi{i}] = dgmgraphs(snaps{i+1},snaps{i},dtarget,...
        slowness,dtrigger,period);
end
[graphsi{numSnaps},stepsi{numSnaps}] = dgmgraphs(snaps{numSnaps},...
    snaps{numSnaps},dtarget,slowness,dtrigger,period);

% Make the step indices cumulative.
for i = 2:numSnaps
    stepsi{i} = stepsi{i} + stepsi{i-1}(end);
end

% Compute clustering coefficients of interpolating graphs.
aci = cell(numSnaps,1);
cci = cell(numSnaps,1);
gci = cell(numSnaps,1);  
for i = 1:numSnaps
    len = length(graphsi{i});
    aci{i} = zeros(1,len);
    cci{i} = zeros(numNodes,len);
    gci{i} = zeros(1,len);
    for j = 1:length(graphsi{i})
        [aci{i}(j),cci{i}(:,j),gci{i}(j)] = avgClusteringCoefficient(graphsi{i}{j});
    end
end

% Draw 4 equally spaced graphs from the second-to-last snapshot to the last
% snapshot.
len = length(graphsi{numSnaps-1});
%xdata = repmat(1:4,1,8);
%ydata = repelem(8:-1:1,4);

figure
p4 = plot(graph(graphsi{end-1}{end}),'LineWidth',2);
xdata = p4.XData;
ydata = p4.YData;
labelnode(p4,1:numNodes,'');
set(gca,'xtick',[],'ytick',[]);
title('Target graph')
set(gca,'fontsize',20)
saveas(gcf,'vdbuntgraph4.eps','epsc')
saveas(gcf,'vdbuntgraph4.fig')

figure
p1 = plot(graph(graphsi{end-1}{1}),'XData',xdata,'YData',ydata,'LineWidth',2);
labelnode(p1,1:numNodes,'');
set(gca,'xtick',[],'ytick',[]);
title('Starting graph')
set(gca,'fontsize',20)
saveas(gcf,'vdbuntgraph1.eps','epsc')
saveas(gcf,'vdbuntgraph1.fig')

figure
p2 = plot(graph(graphsi{end-1}{1+floor(len/3)}),'XData',xdata,'YData',ydata,'LineWidth',2);
labelnode(p2,1:numNodes,'');
set(gca,'xtick',[],'ytick',[]);
title(['Step ' num2str(floor(len/3))])
set(gca,'fontsize',20)
saveas(gcf,'vdbuntgraph2.eps','epsc')
saveas(gcf,'vdbuntgraph2.fig')

figure
p3 = plot(graph(graphsi{end-1}{1+2*floor(len/3)}),'XData',xdata,'YData',ydata,'LineWidth',2);
labelnode(p3,1:numNodes,'');
set(gca,'xtick',[],'ytick',[]);
title(['Step ' num2str(2*floor(len/3))])
set(gca,'fontsize',20)
saveas(gcf,'vdbuntgraph3.eps','epsc')
saveas(gcf,'vdbuntgraph3.fig')

% EXTRAPOLATE
% Uniform attachment/detachment
acu = cell(numSnaps,1);  % clustering coefficients of extrapolating graphs
ccu = cell(numSnaps,1);
gcu = cell(numSnaps,1);
stepsu = cell(numSnaps,1);  % step indices
for i = 1:numSnaps-1
    [acu{i},ccu{i},gcu{i},stepsu{i}] = extrapunif(snaps{i+1},snaps{i},period);
    % Kludge to not connect the graph between snapshots.
    acu{i}(end+1) = NaN;
    gcu{i}(end+1) = NaN;
    stepsu{i}(end+1) = stepsu{i}(end);
end
[acu{numSnaps},ccu{numSnaps},gcu{numSnaps},stepsu{numSnaps}] = extrapunif(snaps{numSnaps},...
    snaps{numSnaps},period);
% Make the starts of the step indices align with the interpolation.
for i = 2:numSnaps
    stepsu{i} = stepsu{i} + stepsi{i}(1);
end

% Preferential attachment/detachment
acp = cell(numSnaps,1);  % clustering coefficients of extrapolating graphs
ccp = cell(numSnaps,1);
gcp = cell(numSnaps,1);
stepsp = cell(numSnaps,1);  % step indices
for i = 1:numSnaps-1
    [acp{i},ccp{i},gcp{i},stepsp{i}] = extrappref(snaps{i+1},snaps{i},period);
    % Kludge to not connect the graph between snapshots.
    acp{i}(end+1) = NaN;
    gcp{i}(end+1) = NaN;
    stepsp{i}(end+1) = stepsp{i}(end);
end
[acp{numSnaps},ccp{numSnaps},gcp{numSnaps},stepsp{numSnaps}] = extrappref(snaps{numSnaps},...
    snaps{numSnaps},period);
% Make the starts of the step indices align with the interpolation.
for i = 2:numSnaps
    stepsp{i} = stepsp{i} + stepsi{i}(1);
end

% Triangle closing
act = cell(numSnaps,1);  % clustering coefficients of extrapolating graphs
cct = cell(numSnaps,1);
gct = cell(numSnaps,1);
stepst = cell(numSnaps,1);  % step indices
for i = 1:numSnaps-1
    [act{i},cct{i},gct{i},stepst{i}] = extraptri(snaps{i+1},snaps{i},period);
    % Kludge to not connect the graph between snapshots.
    act{i}(end+1) = NaN;
    gct{i}(end+1) = NaN;
    stepst{i}(end+1) = stepst{i}(end);
end
[act{numSnaps},cct{numSnaps},gct{numSnaps},stepst{numSnaps}] = extraptri(snaps{numSnaps},...
    snaps{numSnaps},period);
% Make the starts of the step indices align with the interpolation.
for i = 2:numSnaps
    stepst{i} = stepst{i} + stepsi{i}(1);
end

% SAVE
save('vdbuntanalysis.mat','graphsi','stepsi','aci','cci','gci',...
    'stepsu','acu','ccu','gcu',...
    'stepsp','acp','ccp','gcp',...
    'stepst','act','cct','gct',...
    'period')

% PLOT
figure
hold on
p1 = plot(horzcat(stepsi{:}),horzcat(aci{:}),'LineWidth',2);
p2 = plot(horzcat(stepsu{:}),horzcat(acu{:}),'--','LineWidth',1);
p3 = plot(horzcat(stepsp{:}),horzcat(acp{:}),'--','LineWidth',1);
p4 = plot(horzcat(stepst{:}),horzcat(act{:}),'--','LineWidth',1);
scatter(0,0,'r')
for i = 2:numSnaps
    scatter(stepsi{i}(1),aci{i}(1),'r')
end
legend([p1 p2 p3 p4],{'Interpolation','Extrapolation (uniform)',...
    'Extrapolation (preferential)','Extrapolation (triangle closing)'},...
    'Location','best')
title('University friendships')
xlabel('Steps')
ylabel('Mean clustering coefficient')
set(gca,'fontsize',20)
saveas(gcf,'vdbuntac.eps','epsc')
saveas(gcf,'vdbuntac.fig')

figure
hold on
p1 = plot(horzcat(stepsi{:}),horzcat(gci{:}),'LineWidth',2);
p2 = plot(horzcat(stepsu{:}),horzcat(gcu{:}),'--','LineWidth',1);
p3 = plot(horzcat(stepsp{:}),horzcat(gcp{:}),'--','LineWidth',1);
p4 = plot(horzcat(stepst{:}),horzcat(gct{:}),'--','LineWidth',1);
scatter(0,0,'r')
for i = 2:numSnaps
    scatter(stepsi{i}(1),gci{i}(1),'r')
end
legend([p1 p2 p3 p4],{'Interpolation','Extrapolation (uniform)',...
    'Extrapolation (preferential)','Extrapolation (triangle closing)'},...
    'Location','best')
title('University friendships')
xlabel('Steps')
ylabel('Global clustering coefficient')
set(gca,'fontsize',20)
saveas(gcf,'vdbuntgc.eps','epsc')
saveas(gcf,'vdbuntgc.fig')