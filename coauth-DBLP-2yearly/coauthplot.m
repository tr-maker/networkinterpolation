% Plot the mean and global clustering coefficients of the interpolation and 
% extrapolations for the coauthor dataset.
% Save the figure with the mean clustering coefficients as 'coauthac.fig' 
% and 'coauthac.eps'.
% Save the figure with the global clustering coefficients as 'coauthgc.fig' 
% and 'coauthgc.eps'.
% 
% Dependendencies: coauthsnaps.mat, coauthinterp.mat, coauthextrapunif.mat,
% coauthextrappref.mat, coauthextraptri.mat

load('coauthsnaps.mat','snaps')  % cell array of adjacency matrix snapshots
load('coauthinterp.mat','aci','gci','stepsi')
load('coauthextrapunif.mat','acu','gcu','stepsu')
load('coauthextrappref.mat','acp','gcp','stepsp')
load('coauthextraptri.mat','act','gct','stepst')
numSnaps = length(snaps);
numNodes = size(snaps{1},1);

% PLOT
figure
hold on
p1 = plot(horzcat(stepsi{:}),horzcat(aci{:}),'LineWidth',1.25);
p2 = plot(horzcat(stepsu{:}),horzcat(acu{:}),'--','LineWidth',1);
p3 = plot(horzcat(stepsp{:}),horzcat(acp{:}),'--','LineWidth',1);
p4 = plot(horzcat(stepst{:}),horzcat(act{:}),'--','LineWidth',1);
for i = 1:numSnaps
    scatter(stepsi{i}(1),aci{i}(1),'r')
end
legend([p1 p2 p3 p4],{'Interpolation','Extrapolation (uniform)',...
    'Extrapolation (preferential)','Extrapolation (triangle closing)'},...
    'Location','best')
title('Coauthor network')
xlabel('Steps')
ylabel('Mean clustering coefficient')
set(gca,'fontsize',20)
saveas(gcf,'coauthac.eps','epsc')
saveas(gcf,'coauthac.fig')

figure
hold on
p1 = plot(horzcat(stepsi{:}),horzcat(gci{:}),'LineWidth',1.25);
p2 = plot(horzcat(stepsu{:}),horzcat(gcu{:}),'--','LineWidth',1);
p3 = plot(horzcat(stepsp{:}),horzcat(gcp{:}),'--','LineWidth',1);
p4 = plot(horzcat(stepst{:}),horzcat(gct{:}),'--','LineWidth',1);
for i = 1:numSnaps
    scatter(stepsi{i}(1),gci{i}(1),'r')
end
legend([p1 p2 p3 p4],{'Interpolation','Extrapolation (uniform)',...
    'Extrapolation (preferential)','Extrapolation (triangle closing)'},...
    'Location','best')
title('Coauthor network')
xlabel('Steps')
ylabel('Global clustering coefficient')
set(gca,'fontsize',20)
saveas(gcf,'coauthgc.eps','epsc')
saveas(gcf,'coauthgc.fig')