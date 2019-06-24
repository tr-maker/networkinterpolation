% Run the preferential (Barabasi-Albert) attachment/detachment model.
% Return arrays of the mean and global clustering coefficients of the 
% graphs, and an array of the step indices.
% B is the adjacency matrix of the target graph.
% A is the adjacency matrix of the starting graph.
% period is the sampling period (higher values require less storage).
%
% Note: The first entry in the outputs is the zeroth step of the model,
% and the last entry is the last step of the model.

function [ac,cc,gc,steps] = extrappref(B,A,period)
    numNodes = size(A,1);
    
    g = A;  % the changing graph
    numNewEdges = (nnz(B) - nnz(A))/2; % number of edges to add/delete
    if numNewEdges >= 0
        % Add edges.
        [ac,cc,gc] = avgClusteringCoefficient(g);
        steps = 0;
        for j = 1:numNewEdges
            % Pick a uniformly random node and connect it to a random node weighted by degree.
            d = sum(g);
            v1 = randsample(numNodes,1);
            % Check if there exists a node with nonzero degree to connect to.
            potentialconn1 = ~g(v1,:);
            potentialconn1(v1) = 0;
            w = potentialconn1 .* d;
            if nnz(w) == 0
                w = potentialconn1;
                v2 = randsample(numNodes,1,true,w);
                while v1 == v2
                    v2 = randsample(numNodes,1,true,w);
                end
            else
                v2 = randsample(numNodes,1,true,w);
                while v1 == v2
                    v2 = randsample(numNodes,1,true,w);
                end
            end
            %disp(['v1: ', num2str(v1)])
            %disp(['v2: ', num2str(v2)])
            assert(g(v1,v2) == 0);
            assert(g(v2,v1) == 0);
            g(v1,v2) = 1;
            g(v2,v1) = 1;
            if mod(j, period) == 0 || j == numNewEdges
                [ac(end+1), cc(:,end+1), gc(end+1)] = avgClusteringCoefficient(g);
                steps(end+1) = j;
            end
        end
    else
        % Delete edges.
        numNewEdges = -numNewEdges;
        [ac,cc,gc] = avgClusteringCoefficient(g);
        steps = 0;
        for j = 1:numNewEdges
            % Pick a uniformly random node and delete one of its neighbors weighted by degree.
            v1 = randsample(numNodes,1);
            while isempty(find(g(:,v1),1))
                v1 = randsample(numNodes,1);
            end
            nbr1 = find(g(:,v1));
            dinv = 1./sum(g(:,nbr1));
            ind = randsample(length(nbr1),1,true,dinv);
            v2 = nbr1(ind);
            assert(g(v1,v2) == 1);
            assert(g(v2,v1) == 1);
            g(v1,v2) = 0;
            g(v2,v1) = 0;
            if mod(j, period) == 0 || j == numNewEdges
                [ac(end+1), cc(:,end+1), gc(end+1)] = avgClusteringCoefficient(g);
                steps(end+1) = j;
            end
        end
    end
end