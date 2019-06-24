% Plot quartile box plots of mean and global clustering coefficients vs. 
% edit distance for the interpolation with 2 consecutive snapshots from the 
% van de Bunt social network datafiles. The parameter indexA specifies the 
% index of the starting snapshot.
% Save the box plot with the mean clustering coefficients as 
% 'vdbuntenvelopesnap[i]to[i+1]ac.fig' and
% 'vdbuntenvelopesnap[i]to[i+1]ac.eps'.
% Save the box plot with the global clustering coefficients as 
% 'vdbuntenvelopesnap[i]to[i+1]gc.fig' and 
% 'vdbuntenvelopesnap[i]to[i+1]gc.eps'.
%
% Dependencies: vdbuntsnaps.mat

load('vdbuntsnaps.mat','snaps')  % cell array of adjacency matrix snapshots
%numSnaps = length(snaps);
%numNodes = size(snaps{1},1);

indexA = 4;  % SNAP INDEX OF THE STARTING GRAPH.
A = snaps{indexA};  % the starting graph
B = snaps{indexA+1};  % the ending graph
dtarget = 0;  % target edit distance
slowness = 1;  % rate of approach to the target edit distance (higher values are slower)
dtrigger = 0;  % graph distance to the background graph that stops the dynamic graph model
period = 1;  % sampling period (higher values require less storage).
numTrials = 10000;

editDist = nnz(triu(B-A,1));  % edit distance between A and B
acat = cell(editDist+1,1);  % acat{i+1} is the set of clustering coefficients at edit distance i
gcat = cell(editDist+1,1);
for i = 1:numTrials
    [ac,cc,gc,ds] = dgmenvelope(B,A,dtarget,slowness,dtrigger,period);
    len = length(ds);
    for j = 1:length(ds)
        acat{ds(j)+1}(end+1) = ac(j);
        gcat{ds(j)+1}(end+1) = gc(j);
    end
end
quartilesacat = zeros(3,editDist+1); % the ith row is the ith quartile at edit distance i
for i = 1:editDist+1
    quartilesacat(:,i) = quantile(acat{i},[0.25, 0.5, 0.75])';
end
quartilesgcat = zeros(3,editDist+1); % the ith row is the ith quartile at edit distance i
for i = 1:editDist+1
    quartilesgcat(:,i) = quantile(gcat{i},[0.25, 0.5, 0.75])';
end

% PLOT
figure
hold on
x = 0:editDist;
y1 = quartilesacat(1,:);
y2 = quartilesacat(2,:);
y3 = quartilesacat(3,:);
%legend([plot(x,y1),plot(x,y2),plot(x,y3)],{'1st quartile','Median','3rd quartile'},'AutoUpdate','off')
lightblue = [0, 0.4470, 0.7410] * 255;
tintfactor = 0.5;
lightblue = (lightblue + (255 - lightblue) * tintfactor)/255;
fill([x fliplr(x)],[y1 fliplr(y2)],lightblue)
fill([x fliplr(x)],[y2 fliplr(y3)],lightblue)
plot(x,y2,'Color','r','LineWidth',2)
xlim([0 editDist])
ylim([0 0.8])
%title('Quartile box plot of clustering coefficients')
xlabel('Edit distance to target')
ylabel('Mean clustering coefficient')
set(gca,'fontsize',20)
saveas(gcf,['vdbuntenvelopesnap' num2str(indexA) 'to' num2str(indexA+1) 'ac.eps'],'epsc')
saveas(gcf,['vdbuntenvelopesnap' num2str(indexA) 'to' num2str(indexA+1) 'ac.fig'])

figure
hold on
x = 0:editDist;
y1 = quartilesgcat(1,:);
y2 = quartilesgcat(2,:);
y3 = quartilesgcat(3,:);
%legend([plot(x,y1),plot(x,y2),plot(x,y3)],{'1st quartile','Median','3rd quartile'},'AutoUpdate','off')
lightblue = [0, 0.4470, 0.7410] * 255;
tintfactor = 0.5;
lightblue = (lightblue + (255 - lightblue) * tintfactor)/255;
fill([x fliplr(x)],[y1 fliplr(y2)],lightblue)
fill([x fliplr(x)],[y2 fliplr(y3)],lightblue)
plot(x,y2,'Color','r','LineWidth',2)
xlim([0 editDist])
ylim([0 0.8])
%title('Quartile box plot of clustering coefficients')
xlabel('Edit distance to target')
ylabel('Global clustering coefficient')
set(gca,'fontsize',20)
saveas(gcf,['vdbuntenvelopesnap' num2str(indexA) 'to' num2str(indexA+1) 'gc.eps'],'epsc')
saveas(gcf,['vdbuntenvelopesnap' num2str(indexA) 'to' num2str(indexA+1) 'gc.fig'])

% Helper function for creating envelopes.
% Run our graph interpolation until the specified graph distance to the background graph is reached.
% Return the clustering coefficients of the graphs and an array of the edit distances.
% B is the adjacency matrix of the background graph.
% A is the adjacency matrix of the initial current graph.
% dtarget is the target edit distance between the background graph and the current graph.
% slowness is the rate of approach to the target edit distance (higher values are slower).
% dtrigger is the graph distance to the background graph that stops the dynamic graph model.
% period is the sampling period (higher values require less storage).
% 
% Notes: The first entry in the outputs is the zeroth step of the model, 
% and the last entry is the last step of the model.
% False edges are allowed.
% The implementation is efficient provided that A and B are sparse.
% It assumes A and B are undirected.

function [ac,cc,gc,ds] = dgmenvelope(B,A,dtarget,slowness,dtrigger,period)
    n = size(A,1);
    
    % Build upper triangular adjacency matrix of advancing edges.
    % An entry of 1 indicates that adding the corresponding edge is an advancing move.
    % An entry of -1 indicates that deleting the corresponding edge is an advancing move.
    % An entry of 0 indicates that the corresponding edge is not advancing.
    U = triu(B-A,1);
    
    % upper adjacency matrix of current graph
    graph = triu(A,1);
    % current number of advancing moves (edit distance between background graph and current graph)
    d = nnz(U);
    % current number of regressing moves
    duseless = n*(n-1)/2 - d;
    
    %%assert(d == sum(graph ~= triu(B,1),'all'));
    
    % the output
    [ac,cc,gc] = avgClusteringCoefficient(graph + graph');
    ds = d;
    
    if d <= dtrigger
        return;
    end
    
    i = 1;
    while true
        % Make a random move according to the ``Markov'' chain. 
        prb = NaN;
        if d == 0
            prb = 0;
        elseif duseless == 0
            prb = 1;
        else
            % Use a logistic function for d -> p.
            prb = 1/(1 + exp((-d + dtarget)/slowness));  
        end
        
        if binornd(1,prb)
            % Do an advancing move.
            usefulIndices = find(U);
            usefulEdge = usefulIndices(randi(d));
            if U(usefulEdge) == 1
                % Add the edge.
                %%assert(graph(usefulEdge) == 0);
                graph(usefulEdge) = 1;
                U(usefulEdge) = 0;
                d = d - 1;
                duseless = duseless + 1;
            else
                %%assert(U(usefulEdge) == -1)
                % Delete the edge.
                %%assert(graph(usefulEdge) == 1);
                graph(usefulEdge) = 0;
                U(usefulEdge) = 0;
                d = d - 1;
                % The following happens because we allow false edges.
                duseless = duseless + 1;
            end
        else 
            % Do a regressing move (via rejection sampling).
            uselessEdge = sort(randperm(n,2));
            while U(uselessEdge(1),uselessEdge(2)) ~= 0
                uselessEdge = sort(randperm(n,2));
            end
            
            if B(uselessEdge(1),uselessEdge(2)) == 1
                % Delete the edge.
                %%assert(graph(uselessEdge(1),uselessEdge(2)) == 1);
                graph(uselessEdge(1),uselessEdge(2)) = 0;
                U(uselessEdge(1),uselessEdge(2)) = 1;
                d = d + 1;
                duseless = duseless - 1;
            else
                % Add the edge. 
                % The following happens because we allow false edges.
                %%assert(graph(uselessEdge(1),uselessEdge(2)) == 0);
                graph(uselessEdge(1),uselessEdge(2)) = 1;
                U(uselessEdge(1),uselessEdge(2)) = -1;
                d = d + 1;
                duseless = duseless - 1;
            end
        end
        
        %%disp(usefuledges)
        %%disp(uselessedges)
        %%disp(graph)
        %%disp(i)
        %%assert(d == sum(graph ~= triu(B,1),'all'))
        
        if mod(i, period) == 0 || d <= dtrigger
            %%disp(i)
            [ac(end+1), cc(:,end+1), gc(end+1)] = avgClusteringCoefficient(graph + graph');
            ds(end+1) = d;
            if d <= dtrigger
                return;
            end
        end
        i = i + 1;
    end
end