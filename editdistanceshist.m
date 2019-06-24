% Using our interpolation model, plot a histogram of distribution of edit 
% distances over time and compare with the analytic approximation.
% The starting graph is an Erdos-Renyi random graph and the background
% graph is a 2-block SBM graph.
% Do the histogram for two different rate parameters (set to 1 and 10,
% respectively).
% Save the first histogram as 'editdistanceshists1.fig' and
% 'editdistanceshists1.eps'.
% Save the second histogram as 'editdistanceshists10.fig' and
% 'editdistanceshists10.eps'.
% 
% Dependencies: sbm.m, dgmrun.m

% Parameters.
n = 50;  % total number of nodes; should be divisible by 2 (eg 50)
pr = 0.5; % Erdos-Renyi edge connection probability (eg 0.5)
p = 0.9;  % SBM in-cluster connection probability (eg 0.9)
q = 0.1;  % SBM out-of-cluster connection probability (eg 0.1)

dtarget = 10;  % target edit distance
numsteps = 20000;  % number of steps performed
period = 1;  % sampling period (higher values require less storage)
stepsbegin = 10001;  % number of steps where histogram begins

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

% RATE 1 %
slowness = 1;  % rate of approach to the target edit distance (higher values are slower)
[ds,~] = dgmrun(B,A,dtarget,slowness,numsteps,period);
ds = ds(ceil(stepsbegin/period) + 1:end);
figure
histogram(ds,'Normalization','probability');
hold on

% Analytic approximation:
v = zeros(1,2*dtarget+1);  % v(i+1) is the frequency at edit distance i
i = 0:100;
denom = 4*sum(exp(-i.*(i+1)/(2*slowness)));
for k = 0:dtarget
    v(dtarget + 1 + k) = ...
        (exp(-k*(k-1)/(2*slowness)) + exp(-k*(k+1)/(2*slowness)))/denom;
    v(dtarget + 1 - k) = v(dtarget + 1 + k);
end
%%disp(sum(v))  % should be pretty much 1
%%disp(v(1))  % should be pretty much 0
%%disp(v(end))  % should be pretty much 0
plot(0:2*dtarget, v, '-x', 'LineWidth', 2)

title(['Time-averaged limiting distribution' newline 'of edit distances'])
xlabel('Edit distance')
ylabel('Frequency')
legend('Experiment',['Analytic' newline 'approximation']);
set(gca,'fontsize', 20)
saveas(gcf,['editdistanceshists' num2str(slowness) '.eps'],'epsc')
saveas(gcf,['editdistanceshists' num2str(slowness) '.fig'])

% RATE 2 %
slowness = 10;
[ds,~] = dgmrun(B,A,dtarget,slowness,numsteps,period);
ds = ds(ceil(stepsbegin/period) + 1:end);
figure
histogram(ds,'Normalization','probability');
hold on

% Analytic approximation:
v = zeros(1,2*dtarget+1);  % v(i+1) is the frequency at edit distance i
i = 0:100;
denom = 4*sum(exp(-i.*(i+1)/(2*slowness)));
for k = 0:dtarget
    v(dtarget + 1 + k) = ...
        (exp(-k*(k-1)/(2*slowness)) + exp(-k*(k+1)/(2*slowness)))/denom;
    v(dtarget + 1 - k) = v(dtarget + 1 + k);
end
%%disp(sum(v))  % should be pretty much 1
%%disp(v(1))  % should be pretty much 0
%%disp(v(end))  % should be pretty much 0
plot(0:2*dtarget, v, '-x', 'LineWidth', 2)
ylim([0 0.145])

title(['Time-averaged limiting distribution' newline 'of edit distances'])
xlabel('Edit distance')
ylabel('Frequency')
legend('Experiment',['Analytic' newline 'approximation']);
set(gca,'fontsize', 20)
saveas(gcf,['editdistanceshists' num2str(slowness) '.eps'],'epsc')
saveas(gcf,['editdistanceshists' num2str(slowness) '.fig'])