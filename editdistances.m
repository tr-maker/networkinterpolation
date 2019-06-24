% Using our interpolation model, plot the edit distance over time for 
% different rate parameters. 
% The starting graph is an Erdos-Renyi random graph and the background
% graph is a 2-block SBM graph.
% Save the plot as 'editdistances.fig' and 'editdistances.eps'.
%
% Dependencies: sbm.m, dgmrun.m

% Parameters
n = 50;  % total number of nodes; should be divisible by 2 (eg 50)
pr = 0.5; % Erdos-Renyi edge connection probability (eg 0.5)
p = 0.9;  % SBM in-cluster connection probability (eg 0.9)
q = 0.1;  % SBM out-of-cluster connection probability (eg 0.1)

dtarget = 10;  % target edit distance
slownesses = [1, 20, 100];  % rates of approach to the target edit distance (higher values are slower)
numsteps = 2000;  % number of steps performed
period = 1;  % sampling period (higher values require less storage)

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

% starting graph
A = sbm(repelem(1, n), pr, pr);
% background graph
B = sbm(repelem(1:2, n/2), p, q);
figure
for i = 1:length(slownesses)
    hold on
    [ds,~] = dgmrun(B,A,dtarget,slownesses(i),numsteps,period);
    plot(0:period:numsteps, ds, 'LineWidth', 2)
end

title(['Evolution of graph edit distance (d_t = ' num2str(dtarget) ')'])
xlabel('Steps')
ylabel('Graph edit distance')
leg = cell(1,length(slownesses));
for i = 1:length(slownesses)
    leg{i} = ['s = ' num2str(slownesses(i))];
end
legend(leg);
set(gca,'fontsize', 20)
saveas(gcf,'editdistances.fig')
saveas(gcf,'editdistances.eps','epsc')