% Consider two scenarios where a 2-block graph evolves into a 3-block
% graph:
% 1. The second block in the 2-block graph splits into 2 equally sized 
% blocks.
% 2. The blocks in the 3-block graph have no a priori relation to the 2- 
% block graph.
% Use our model to do interpolation in both scenarios.
% Save and plot various statistics:
% 1. Recovery rate using our QR algorithm and the singular value distance 
% to the target graph. The singular value distance is between A and B is 
% sqrt(1 - sigma_min(A^TB)^2) if A and B have orthonormal columns.
% Plot the recovery rate and singular value distance on the same figure.
% Save the figure as 'twotothreeblockrecovery1.fig' and
% 'twotothreeblockrecovery1.eps' for Scenario 1, and 
% 'twotothreeblockrecovery2.fig' and 'twotothreeblockrecovery2.eps' for
% Scenario 2.
% 2. The spectrum of the symmetrically normalized adjacency matrices of 
% the interpolating graphs.
% Save the figure as 'twotothreeblockeigenvalues1.fig' and 
% 'twotothreeblockeigenvalues1.eps' for Scenario 1, and 
% 'twotothreeblockeigenvalues2.fig' and 'twotothreeblockeigenvalues2.eps' 
% for Scenario 2.
% 3. The spectrum of the equispaced linear interpolation between the 
% symmetrically normalized adjacency matrices of the starting and target 
% graphs. 
% Save the figure as 'twotothreeblockeigenvaluesalg1.fig' and 
% 'twotothreeblockeigenvaluesalg1.eps' for Scenario 1, and 
% 'twotothreeblockeigenvaluesalg2.fig' and 
% 'twotothreeblockeigenvaluesalg2.eps' for Scenario 2.
%
% Save all data relating to 1 and 2 above as 'twotothreeblock.mat'.
% Save all data relating to 3 above as 'twotothreeblockcomparison.mat'.
% 
% Dependencies: dgmtrigger.m, sbm.m, numcorrect.m

% Parameters.
n = 120;  % total number of nodes; should be divisible by 4 and 3 (eg 120)
p = 0.9;  % in-cluster connection probability (eg 0.9)
q = 0.1;  % out-of-cluster connection probability (eg 0.1)

dtarget = 0;  % target edit distance
slowness = 1;  % rate of approach to the target edit distance (higher values are slower)
dtrigger = 0;  % graph distance to the target graph that stops the dynamic graph model
period = 1;  % sampling period (higher values require less storage)

% Check parameters.
if mod(n,12) ~= 0
    error('Total number of nodes, n, is not divisible by 12.')
end
alpha = p*n/log(n);
beta = q*n/log(n);
if sqrt(alpha) - sqrt(beta) <= sqrt(3)
    error(['Connection probabilities are not in the theoretical exact '...
        'recovery regime for 3 clusters.'])
end

% starting graph
truthA = repelem(1:2, n/2);
A = sbm(truthA, p, q);
% target graph in scenario 1
truthB1 = [repelem(1, n/2), repelem(2:3, n/4)];
B1 = sbm(truthB1, p, q);
% target graph in scenario 2
truthB2 = repelem(1:3, n/3);
truthB2 = truthB2(randperm(n));
B2 = sbm(truthB2, p, q);

[graphs1,numSteps1] = dgmtrigger(B1,A,dtarget,slowness,dtrigger,period);
[graphs2,numSteps2] = dgmtrigger(B2,A,dtarget,slowness,dtrigger,period);
assert(isequal(B1, graphs1{end}));
assert(isequal(B2, graphs2{end}));
len1 = length(graphs1);
len2 = length(graphs2);
% Normalize the adjacency matrices.
for i = 1:len1
    dtotheneghalf = 1./sqrt(sum(graphs1{i}));
    dtotheneghalf(~isfinite(dtotheneghalf)) = 0;
    dtotheneghalf = diag(dtotheneghalf);
    graphs1{i} = dtotheneghalf * graphs1{i} * dtotheneghalf;
end
for i = 1:len2
    dtotheneghalf = 1./sqrt(sum(graphs2{i}));
    dtotheneghalf(~isfinite(dtotheneghalf)) = 0;
    dtotheneghalf = diag(dtotheneghalf);
    graphs2{i} = dtotheneghalf * graphs2{i} * dtotheneghalf;
end
[v1,~] = eigs(graphs1{end},3);
[v2,~] = eigs(graphs2{end},3);

recovery1 = zeros(1,len1); % QR recovery rate
svdist1 = zeros(1,len1); % singular value distance to target graph
eigenvalues1 = zeros(n,len1); % eigenvalues
recovery2 = zeros(1,len2);
svdist2 = zeros(1,len2);
eigenvalues2 = zeros(n,len2);
for i = 1:len1
    [v,~] = eigs(graphs1{i}, 3);
    % piv consists of the first k entries in the permutation vector 
    % determined by CPQR
    [~,~,piv] = qr(v','vector');
    piv = piv(1:3);
    % u is the orthogonal matrix in the polar decomposition of v'
    [ut,~,vt] = svd(v(piv,:)');
    u = ut * vt';
    % return the index of the maximum value in each column of abs(u' * v')
    [~,c] = max(abs(u' * v'));
    
    recovery1(i) = numcorrect(c, truthB1)/n;
    svdist1(i) = sqrt(1 - svds(v'*v1,1,'smallest')^2);
    eigenvalues1(:,i) = sort(eig(graphs1{i}),'descend');
end
for i = 1:len2
    [v,~] = eigs(graphs2{i}, 3);
    % piv consists of the first k entries in the permutation vector 
    % determined by CPQR
    [~,~,piv] = qr(v','vector');
    piv = piv(1:3);
    % u is the orthogonal matrix in the polar decomposition of v'
    [ut,~,vt] = svd(v(piv,:)');
    u = ut * vt';
    % return the index of the maximum value in each column of abs(u' * v')
    [~,c] = max(abs(u' * v'));
    
    recovery2(i) = numcorrect(c, truthB2)/n;
    svdist2(i) = sqrt(1 - svds(v'*v2,1,'smallest')^2);
    eigenvalues2(:,i) = sort(eig(graphs2{i}),'descend');
end

save('twotothreeblock.mat','n','p','q',...
    'dtarget','slowness','dtrigger','period',...
    'graphs1','graphs2','recovery1','svdist1','eigenvalues1',...
    'recovery2','svdist2','eigenvalues2')

figure
plot([0:period:period*(len1-2), numSteps1], recovery1, 'LineWidth', 2)
hold on
plot([0:period:period*(len1-2), numSteps1], svdist1, 'LineWidth', 2)
%title('2-to-3 block graph, Scenario 1')
xlabel('Steps')
xticks([0 500 1000 1500 2000])
yticks([0 0.2 0.4 0.6 0.8 1])
legend('Recovery rate','Subspace distance','Location','southwest')
set(gca,'fontsize',20)
saveas(gcf,'twotothreeblockrecovery1.eps','epsc')
saveas(gcf,'twotothreeblockrecovery1.fig')

figure
plot([0:period:period*(len2-2), numSteps2], recovery2, 'LineWidth', 2)
hold on
plot([0:period:period*(len2-2), numSteps2], svdist2, 'LineWidth', 2)
%title('2-to-3 block graph, Scenario 2')
xlabel('Steps')
xticks([0 1000 2000 3000 4000])
yticks([0 0.2 0.4 0.6 0.8 1])
legend('Recovery rate','Subspace distance','Location','southwest')
set(gca,'fontsize',20)
saveas(gcf,'twotothreeblockrecovery2.eps','epsc')
saveas(gcf,'twotothreeblockrecovery2.fig')

figure
hold on
x = [0:period:period*(len1-2), numSteps1];
for i = 1:3
    p = plot(x,eigenvalues1(i,:),'LineWidth',2);
end
for i = 4:n
    p = plot(x,eigenvalues1(i,:),'LineWidth',2);
    p.Color(4) = 0.1;  % set alpha value
end
%plot([0:period:period*(len1-2), numSteps1], eigenvalues1', 'LineWidth', 3)
%title('2-to-3 block graph, Scenario 1')
xlabel('Steps')
ylabel('Eigenvalues')
ylim([-0.2 1])
set(gca,'fontsize',20)
saveas(gcf,'twotothreeblockeigenvalues1.eps','epsc')
saveas(gcf,'twotothreeblockeigenvalues1.fig')
figure
hold on
x = [0:period:period*(len2-2), numSteps2];
for i = 1:4
    p = plot(x,eigenvalues2(i,:),'LineWidth',2);
end
for i = 5:n
    p = plot(x,eigenvalues2(i,:),'LineWidth',2);
    p.Color(4) = 0.1;
end
%plot([0:period:period*(len2-2), numSteps2], eigenvalues2', 'LineWidth', 3)
%title('2-to-3 block graph, Scenario 2')
xlabel('Steps')
ylabel('Eigenvalues')
ylim([-0.2 1])
set(gca,'fontsize',20)
saveas(gcf,'twotothreeblockeigenvalues2.eps','epsc')
saveas(gcf,'twotothreeblockeigenvalues2.fig')

% Plot the eigenvalues of the algebraic interpolation.
A = graphs1{1};
B1 = graphs1{end};
B2 = graphs2{end};
eigenvaluescmp1 = zeros(n,len1);
eigenvaluescmp2 = zeros(n,len2);
for i = 0:len1-2
    eigenvaluescmp1(:,i+1) = sort(eig((1 - i*period/numSteps1)*A + (i*period/numSteps1)*B1),'descend');
end
eigenvaluescmp1(:,len1) = sort(eig(B1),'descend');
for i = 0:len2-2
    eigenvaluescmp2(:,i+1) = sort(eig((1 - i*period/numSteps2)*A + (i*period/numSteps2)*B2),'descend');
end
eigenvaluescmp2(:,len2) = sort(eig(B2),'descend');
save('twotothreeblockcomparison.mat','n','p','q','period',...
    'A','B1','B2','eigenvaluescmp1','eigenvaluescmp2')
figure
hold on
x = [0:period:period*(len1-2), numSteps1];
for i = 1:3
    p = plot(x,eigenvaluescmp1(i,:),'LineWidth',2);
end
for i = 4:n
    p = plot(x,eigenvaluescmp1(i,:),'LineWidth',2);
    p.Color(4) = 0.1;
end
%plot([0:period:period*(len1-2), numSteps1], eigenvaluescmp1', 'LineWidth', 3)
%title('2-to-3 block graph, Scenario 1 (algebraic interpolation)')
xlabel('Steps')
ylabel('Eigenvalues')
ylim([-0.2 1])
set(gca,'fontsize',20)
saveas(gcf,'twotothreeblockeigenvaluesalg1.eps','epsc')
saveas(gcf,'twotothreeblockeigenvaluesalg1.fig')
figure
hold on
x = [0:period:period*(len2-2), numSteps2];
for i = 1:4
    p = plot(x,eigenvaluescmp2(i,:),'LineWidth',2);
end
for i = 5:n
    p = plot(x,eigenvaluescmp2(i,:),'LineWidth',2);
    p.Color(4) = 0.1;
end
%plot([0:period:period*(len2-2), numSteps2], eigenvaluescmp2', 'LineWidth', 3)
%title('2-to-3 block graph, Scenario 2 (algebraic interpolation)')
xlabel('Steps')
ylabel('Eigenvalues')
ylim([-0.2 1])
set(gca,'fontsize',20)
saveas(gcf,'twotothreeblockeigenvaluesalg2.eps','epsc')
saveas(gcf,'twotothreeblockeigenvaluesalg2.fig')