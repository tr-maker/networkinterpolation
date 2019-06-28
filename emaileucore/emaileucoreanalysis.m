% Interpolate between the 2 snapshots from the emaileucore datafile.
% Extrapolate to the 2nd snapshot from the emaileucore datafile
% using the uniform attachment model, preferential attachment model, and 
% the triangle closing model.
% Plot mean and global clustering coefficients of real data, interpolation, 
% and extrapolations.
% Save the figure with the mean clustering coefficients as 
% 'emaileucoreac.fig' and 'emaileucoreac.eps'.
% Save the figure with the global clustering coefficients as 
% 'emaileucoregc.fig' and 'emaileucoregc.eps'.
% Save the experiment data in 'emaileucoreanalysis.mat.'
% 
% Dependencies: emaileucoresnaps.mat, dgm.m, avgClusteringCoefficient.m

% cell array of adjacency matrix snapshots
load('emaileucoresnaps.mat','snaps')
numNodes = size(snaps{2},1);

% Interpolation
dtarget = 0;  % target edit distance
slowness = 1;  % rate of approach to the target edit distance (higher values are slower)
dtrigger = 0;  % graph distance to the target graph that stops the dynamic graph model
period = 50;  % sampling period (higher values require less storage)
[aci,cci,gci,stepsi] = dgm(snaps{2},snaps{1},dtarget,...
        slowness,dtrigger,period);
disp('Interpolation done')

% Extrapolation
% UNIFORM ATTACHMENT MODEL
m = 16; % number of nodes that each new node connects to
A = snaps{1};

% Start with a clique on m vertices.
A(1:m,1:m) = ~eye(m);
acu = zeros(1,ceil((numNodes-m)/period)+1);
ccu = zeros(numNodes,ceil((numNodes-m)/period)+1);
gcu = zeros(1,ceil((numNodes-m)/period)+1);
[acu(1),ccu(:,1),gcu(1)] = avgClusteringCoefficient(A);

for i = 1:numNodes-m
    % Each new node connects to m uniformly random nodes.
    randnodes = randsample(1:i+m-1,m);
    A(i+m,randnodes) = 1;
    A(randnodes,i+m) = 1;
    if mod(i,period) == 0
        %disp(i)
        [acu(i/period+1),ccu(:,i/period+1),gcu(i/period+1)] = avgClusteringCoefficient(A);
    end
end
if mod(numNodes-m,period) ~= 0
    [acu(end),ccu(:,end),gcu(end)] = avgClusteringCoefficient(A);
end

xu = 0:period:numNodes-m;
if mod(numNodes-m,period) ~= 0
    xu = [xu numNodes-m];
end
%disp(['Uniform: ' num2str(nnz(A)/2) ' edges'])

% PREFERENTIAL ATTACHMENT MODEL
A = snaps{1};
m = 16; % number of nodes that each new node connects to

% Start with a clique on m vertices.
A(1:m,1:m) = ~eye(m);
acp = zeros(1,ceil((numNodes-m)/period)+1);
ccp = zeros(numNodes,ceil((numNodes-m)/period)+1);
gcp = zeros(1,ceil((numNodes-m)/period)+1);
[acp(1),ccp(:,1),gcp(1)] = avgClusteringCoefficient(A);

for i = 1:numNodes-m
    % Each new node connects to (approximately) m random nodes weighted by degree.
    w = full(sum(A(1:i+m-1,1:i+m-1)));
    randnodes = randsample(1:i+m-1,m,true,w);
    A(i+m,randnodes) = 1;
    A(randnodes,i+m) = 1;
    if mod(i,period) == 0
        %disp(i)
        [acp(i/period+1),ccp(:,i/period+1),gcp(i/period+1)] = avgClusteringCoefficient(A);
    end
end
if mod(numNodes-m,period) ~= 0
    [acp(end),ccp(:,end),gcp(end)] = avgClusteringCoefficient(A);
end

xp = 0:period:numNodes-m;
if mod(numNodes-m,period) ~= 0
    xp = [xp numNodes-m];
end
%disp(['Preferential: ' num2str(nnz(A)/2) ' edges'])

% TRIANGLE CLOSING MODEL
mr = 6;  % number of nodes chosen as parent nodes
pr = 0.5;  % probability of each connection to parent
mn = 26;  % number of nodes chosen from parents' neighbors
pn = 0.5;  % probability of each connection to parents' neighbor (usually the same as pr)

A = snaps{1};
% Start with a clique on mr+mn+1 vertices.
minit = mr+mn+1;
A(1:minit,1:minit) = ~eye(minit);
act = zeros(1,ceil((numNodes-minit)/period)+1);
cct = zeros(numNodes,ceil((numNodes-minit)/period)+1);
gct = zeros(1,ceil((numNodes-minit)/period)+1);
[act(1),cct(:,1),gct(1)] = avgClusteringCoefficient(A);

for i = 1:numNodes-minit
    % Pick mr nodes uniformly at random (without replacement).
    parents = randsample(1:i+minit-1,mr);
    % Connect to each node independently with probability pr.
    for j = 1:mr
        prob = (rand < pr);
        if prob
            A(i+minit,parents(j)) = 1;
            A(parents(j),i+minit) = 1;
        end
    end
    % Pick mn nodes uniformly at random (without replacement) from the
    % parents' neighborhoods.
    parentsneighbors = any(A(parents,:));
    parentsneighbors(parents) = 0;
    parentsneighbors(i+minit) = 0;
    [~,~,parentsneighbors] = find(parentsneighbors);
    k = min(length(parentsneighbors),mn);
    parentsneighbors = randsample(parentsneighbors,k);
    % Connect to each node independently with probability pn.
    for j = 1:k
        prob = (rand < pn);
        if prob
            A(i+minit,parentsneighbors(j)) = 1;
            A(parentsneighbors(j),i+minit) = 1;
        end
    end
    
    if mod(i,period) == 0
        %disp(i)
        [act(i/period+1),cct(:,i/period+1),gct(i/period+1)] = avgClusteringCoefficient(A);
    end
end
if mod(numNodes-minit,period) ~= 0
    [act(end),cct(:,end),gct(end)] = avgClusteringCoefficient(A);
end

xt = 0:period:numNodes-minit;
if mod(numNodes-minit,period) ~= 0
    xt = [xt numNodes-minit];
end
%disp(['Triangle closing: ' num2str(nnz(A)/2) ' edges'])

% SAVE
save('emaileucoreanalysis.mat','stepsi','aci',...
    'acu','ccu','gcu',...
    'acp','ccp','gcp',...
    'act','cct','gct',...
    'period')

% PLOT
figure
hold on
plot(stepsi,aci,'LineWidth',2)
plot(xu,acu,'--','LineWidth',2)
plot(xp,acp,'--','LineWidth',2)
plot(xt,act,'--','LineWidth',2)
line([0 stepsi(end)],[aci(end),aci(end)],'Color','k','LineWidth',1)
xlim([0 stepsi(end)])
ylim([0 0.5])
title('Institution messages')
xlabel('Steps')
ylabel('Mean clustering coefficient')
legend({'Interpolation (s = 1)',...
    'Extrapolation (uniform)','Extrapolation (preferential)',...
    'Extrapolation (triangle closing)'},...
    'Location','best')
set(gca,'fontsize',20)
saveas(gcf,'emaileucoreac.eps','epsc')
saveas(gcf,'emaileucoreac.fig')

figure
hold on
plot(stepsi,gci,'LineWidth',2)
plot(xu,gcu,'--','LineWidth',2)
plot(xp,gcp,'--','LineWidth',2)
plot(xt,gct,'--','LineWidth',2)
line([0 stepsi(end)],[gci(end),gci(end)],'Color','k','LineWidth',1)
xlim([0 stepsi(end)])
ylim([0 0.5])
title('Institution messages')
xlabel('Steps')
ylabel('Global clustering coefficient')
legend({'Interpolation (s = 1)',...
    'Extrapolation (uniform)','Extrapolation (preferential)',...
    'Extrapolation (triangle closing)'},...
    'Location','best')
set(gca,'fontsize',20)
saveas(gcf,'emaileucoregc.eps','epsc')
saveas(gcf,'emaileucoregc.fig')