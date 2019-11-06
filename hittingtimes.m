% Using our interpolation model, plot a histogram of the distribution of 
% hitting times and compare with the analytic approximation.
% The starting graph is an Erdos-Renyi random graph and the target
% graph is a 2-block SBM graph.
% Do the histogram for two different rate parameters (set to 1 and 10,
% respectively).
% Save the first histogram as 'hittingtimess1.fig' and
% 'hittingtimess1.eps'.
% Save the second histogram as 'hittingtimess10.fig' and
% 'hittingtimess10.eps'.
% Dependencies: sbm.m, dgmtrigger.m

% Parameters.
n = 50;  % total number of nodes; should be divisible by 2 (eg 50)
pr = 0.5; % Erdos-Renyi edge connection probability (eg 0.5)
p = 0.9;  % SBM in-cluster connection probability (eg 0.9)
q = 0.1;  % SBM out-of-cluster connection probability (eg 0.1)

dtarget = 10;  % target edit distance
s = 10;  % rate of approach to the target edit distance (higher values are slower)
dtrigger = dtarget;  % graph distance to the target graph that stops the dynamic graph model
period = 1;  % sampling period (higher values require less storage) - actually doesn't matter
numtrials = 1000;

% Check parameters.
if mod(n,2) ~= 0
    error('Total number of nodes, n, is not divisible by 2.')
end
if pr <= log(n)/n
    warning(['Erdos-Renyi connection probability is not in the '... 
        'theoretical connected graph regime.'])
end
alpha = p*n/log(n);
beta = q*n/log(n);
if sqrt(alpha) - sqrt(beta) <= sqrt(2)
    error(['SBM connection probabilities are not in the theoretical '...
        'exact recovery regime for 2 clusters.'])
end

height = 0.7; % height of vertical lines in plot

% starting graph
A = sbm(repelem(1, n), pr, pr);
% target graph
B = sbm(repelem(1:2, n/2), p, q);

% RATE 1 %
% Hitting times
htimes = zeros(1,numtrials);
for i = 1:numtrials
    %disp(i)
    [~,htimes(i)] = dgmtrigger(B,A,dtarget,s,dtrigger,period);
end
figure
p1 = histogram(htimes,'Normalization','probability');
hold on
%line([mean(hittingtimes) mean(hittingtimes)], [0 height])
disp(mean(htimes))
hold on

% Analytic approximation:
d0 = nnz(triu(B-A,1));
%k = 1:(d0 - dtarget);
term = exp(-1/s) * (1 - exp(-(d0-dtarget)/s))/(1 - exp(-1/s));
sums = term;
i = 1;
while term > eps
    i = i + 1;
    term = exp(-sum(1:i)/s) * (1 - exp(-i*(d0-dtarget)/s)) / (1 - exp(-i/s));
    sums = sums + term;
end
disp(i)
h = d0 - dtarget + 2 * sums;
line([h h],[0 height],'Color','k','LineWidth',1)
disp(h)

% PLOT
title('Hitting times')
xlabel('Hitting time')
ylabel('Frequency')
%xl = xlim;
set(gca,'fontsize',20)
saveas(gcf,['hittingtimess' num2str(s) '.eps'],'epsc')
saveas(gcf,['hittingtimess' num2str(s) '.fig'])

% RATE 2 %
% Hitting times
s = 1;
htimes = zeros(1,numtrials);
for i = 1:numtrials
    %disp(i)
    [~,htimes(i)] = dgmtrigger(B,A,dtarget,s,dtrigger,period);
end
figure
p2 = histogram(htimes,'Normalization','probability');
hold on
%line([mean(hittingtimes) mean(hittingtimes)], [0 height])
disp(mean(htimes))
hold on

% Analytic approximation:
d0 = nnz(triu(B-A,1));
%k = 1:(d0 - dtarget);
term = exp(-1/s) * (1 - exp(-(d0-dtarget)/s))/(1 - exp(-1/s));
sums = term;
i = 1;
while term > eps
    i = i + 1;
    term = exp(-sum(1:i)/s) * (1 - exp(-i*(d0-dtarget)/s)) / (1 - exp(-i/s));
    sums = sums + term;
end
disp(i)
h = d0 - dtarget + 2 * sums;
line([h h],[0 height],'Color','k','LineWidth',1)
disp(h)

% PLOT
title('Hitting times')
xlabel('Hitting time')
ylabel('Frequency')
%xlim([-inf xl(2)])
set(gca,'fontsize',20)
saveas(gcf,['hittingtimess' num2str(s) '.eps'],'epsc')
saveas(gcf,['hittingtimess' num2str(s) '.fig'])

% title('Hitting times')
% xlabel('Hitting time')
% ylabel('Frequency')
% legend([p1,p2],{'s = 1','s = 10'})
% %xlim([580 900])
% set(gca,'fontsize',20)
% saveas(gcf,'hittingtimes.eps','epsc')
% saveas(gcf,'hittingtimes.fig')