% Run our graph interpolation model for a specified number of steps.
% Return a cell array of graphs (in adjacency matrix form) and a vector of 
% edit distances between the background graph and current graph.
% b is the adjacency matrix of the background graph.
% a is the adjacency matrix of the initial current graph.
% dtarget is the target edit distance between the background graph and the current graph.
% slowness is the rate of approach to the target edit distance (higher values are slower).
% numsteps is the number of steps performed.
% period is the sampling period (higher values require less storage).
% normalized is a boolean that determines whether the adjacency matrix is normalized (D^(-1/2)AD^(-1/2)).
%
% Notes: The first entry in the outputs is the zeroth step of the model.
% False edges are allowed.
% The implementation is efficient provided that A and B are sparse.
% It assumes A and B are undirected.

function [graphs,ds] = dgmrun(B,A,dtarget,slowness,numsteps,period)
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
    numstorage = ceil(numsteps/period);
    ds = zeros(1,numstorage);  % output vector of edit distances
    graphs = cell(1,numstorage);  % output vector of graphs
    
    ds(1) = d;
    graphs{1} = graph + graph';
    if numsteps == 0
        return;
    end
    
    for i = 1:numsteps
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
        
        if mod(i, period) == 0
            %%disp(i)
            ds(i/period + 1) = d;
            graphs{i/period + 1} = graph + graph';
        end
    end
end