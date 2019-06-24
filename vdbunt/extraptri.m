% Run the triangle closing (Jackson-Rogers) attachment/detachment model.
% Return arrays of the mean and global clustering coefficients of the 
% graphs, and an array of the step indices.
% B is the adjacency matrix of the target graph.
% A is the adjacency matrix of the starting graph.
% period is the sampling period (higher values require less storage).
%
% Note: The first entry in the outputs is the zeroth step of the model,
% and the last entry is the last step of the model.

function [ac,cc,gc,steps] = extraptri(B,A,period)
    numNodes = size(A,1);
    
    % Parameters
    mr = 1;  % number of nodes chosen as parent nodes
    pr = 1;  % probability of each connection to parent
    mn = 1;  % number of nodes chosen from parents' neighbors
    pn = 1;  % probability of each connection to parents' neighbor (usually the same as pr)
    
    g = A;  % the changing graph
    numNewEdges = (nnz(B) - nnz(A))/2; % number of edges to add/delete
    if numNewEdges >= 0
        % Add edges.
        [ac,cc,gc] = avgClusteringCoefficient(g);
        steps = 0;
        if numNewEdges == 0
            return;
        end
        
        numEdgesAdded = 0;
        while true
            % Pick a random node uniformly at random.
            v1 = randsample(numNodes,1);
            % Pick mr nodes uniformly at random (without replacement).
            parents = randsample([1:v1-1 v1+1:numNodes],mr);
            % Connect to each node independently with probability pr.
            for i = 1:mr
                prob = (rand < pr);
                if prob
                    if g(v1,parents(i)) == 1
                        continue;
                    end
                    assert(g(v1,parents(i)) == 0);
                    assert(g(parents(i),v1) == 0);
                    numEdgesAdded = numEdgesAdded + 1;
                    g(v1,parents(i)) = 1;
                    g(parents(i),v1) = 1;
                    if mod(numEdgesAdded, period) == 0 || numEdgesAdded == numNewEdges
                        [ac(end+1), cc(:,end+1), gc(end+1)] = avgClusteringCoefficient(g);
                        steps(end+1) = numEdgesAdded;
                    end
                    if numEdgesAdded >= numNewEdges
                        return;
                    end
                end
            end
            % Pick mn nodes uniformly at random (without replacement) from
            % the parents' neighborhoods.
            parentsneighbors = any(g(parents,:));
            parentsneighbors(parents) = 0;
            parentsneighbors(v1) = 0;
            [~,~,parentsneighbors] = find(parentsneighbors);
            k = min(length(parentsneighbors),mn);
            parentsneighbors = randsample(parentsneighbors,k);
            % Connect to each node independently with probability pn.
            for i = 1:k
                prob = (rand < pn);
                if prob
                    if g(v1,parentsneighbors(i)) == 1
                        continue;
                    end
                    assert(g(v1,parentsneighbors(i)) == 0);
                    assert(g(parentsneighbors(i),v1) == 0);
                    numEdgesAdded = numEdgesAdded + 1;
                    g(v1,parentsneighbors(i)) = 1;
                    g(parentsneighbors(i),v1) = 1;
                    if mod(numEdgesAdded, period) == 0 || numEdgesAdded == numNewEdges
                        [ac(end+1), cc(:,end+1), gc(end+1)] = avgClusteringCoefficient(g);
                        steps(end+1) = numEdgesAdded;
                    end
                    if numEdgesAdded >= numNewEdges
                        return;
                    end
                end
            end
        end
    else
        % Delete edges.
        numNewEdges = -numNewEdges;
        [ac,cc,gc] = avgClusteringCoefficient(g);
        steps = 0;
        for j = 1:numNewEdges
            % Pick a random node uniformly at random.
            d = sum(g);
            v1 = randsample(numNodes,1);
            while d(v1) == 0
                v1 = randsample(numNodes,1);
            end
            
            % Delete an edge uniformly at random that does not belong to a
            % triangle. If there is no such edge, then weight by the number
            % of triangles the edge participates in.
            % The (i,j) entry of the following matrix counts the number of 
            % triangles that edge (i,j) participates in.
            triadj = g*g.*g;
            triadjinv = 1./triadj(v1,:);
            w = triadjinv .* g(v1,:);
            w(isnan(w)) = 0;
            if any(isinf(w))
                w = isinf(w);
            end
            v2 = randsample(numNodes,1,true,w);
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