% Run the uniform attachment/detachment model.
% Return arrays of the mean and global clustering coefficients of the 
% graphs, and an array of the step indices.
% B is the adjacency matrix of the target graph.
% A is the adjacency matrix of the starting graph.
% period is the sampling period (higher values require less storage).
%
% Note: The first entry in the outputs is the zeroth step of the model, 
% and the last entry is the last step of the model.

function [ac,gc,steps] = extrapunif(B,A,period)
    numNodes = size(A,1);

    g = A;  % the changing graph
    numNewEdges = (nnz(B) - nnz(A))/2; % number of edges to add/delete
    if numNewEdges >= 0
        % Add edges.
        % Inline the clustering coefficient calculation.
        deg = sum(g, 2); %Determine node degrees
        cn = diag(g*triu(g)*g); %Number of triangles for each node
        %The local clustering coefficient of each node
        c = zeros(size(deg));
        c(deg > 1) = 2 * cn(deg > 1) ./ (deg(deg > 1).*(deg(deg > 1) - 1)); 
        %Average clustering coefficient of the graph
        %acc = mean(c(deg > 1)); 
        %acc = mean(c(deg > 1),'omitnan');
        ac = mean(c);
        dummy = sum(cn) / sum(triu(g*g,1),'all');
        dummy(isnan(dummy)) = 0;
        gc = dummy;
            
        steps = 0;
        for j = 1:numNewEdges
            % Pick a uniformly random node and connect it to a uniformly random node.
            v1 = randsample(numNodes,1);
            v2 = randsample(numNodes,1);
            while g(v1,v2) == 1 || v1 == v2
                v2 = randsample(numNodes,1);
            end
            assert(g(v1,v2) == 0);
            assert(g(v2,v1) == 0);
            g(v1,v2) = 1;
            g(v2,v1) = 1;
            if mod(j, period) == 0 || j == numNewEdges
                % Inline the clustering coefficient calculation.
                deg = sum(g, 2); %Determine node degrees
                cn = diag(g*triu(g)*g); %Number of triangles for each node
                %The local clustering coefficient of each node
                c = zeros(size(deg));
                c(deg > 1) = 2 * cn(deg > 1) ./ (deg(deg > 1).*(deg(deg > 1) - 1)); 
                %Average clustering coefficient of the graph
                %acc = mean(c(deg > 1)); 
                %acc = mean(c(deg > 1),'omitnan');
                ac(end+1) = mean(c);
                dummy = sum(cn) / sum(triu(g*g,1),'all');
                dummy(isnan(dummy)) = 0;
                gc(end+1) = dummy;
                
                steps(end+1) = j;
                
                disp(j)
            end
        end
    else
        % Delete edges.
        numNewEdges = -numNewEdges;
        % Inline the clustering coefficient calculation.
        deg = sum(g, 2); %Determine node degrees
        cn = diag(g*triu(g)*g); %Number of triangles for each node
        %The local clustering coefficient of each node
        c = zeros(size(deg));
        c(deg > 1) = 2 * cn(deg > 1) ./ (deg(deg > 1).*(deg(deg > 1) - 1)); 
        %Average clustering coefficient of the graph
        %acc = mean(c(deg > 1)); 
        %acc = mean(c(deg > 1),'omitnan');
        ac = mean(c);
        dummy = sum(cn) / sum(triu(g*g,1),'all');
        dummy(isnan(dummy)) = 0;
        gc = dummy;
        
        steps = 0;
        for j = 1:numNewEdges
            % Pick a uniformly random node and delete one of its neighbors uniformly at random.
            v1 = randsample(numNodes,1);
            while isempty(find(g(:,v1),1))
                v1 = randsample(numNodes,1);
            end
            nbr1 = find(g(:,v1));
            ind = randsample(length(nbr1),1);
            v2 = nbr1(ind);
            assert(g(v1,v2) == 1);
            assert(g(v2,v1) == 1);
            g(v1,v2) = 0;
            g(v2,v1) = 0;
            if mod(j, period) == 0 || j == numNewEdges
                % Inline the clustering coefficient calculation.
                deg = sum(g, 2); %Determine node degrees
                cn = diag(g*triu(g)*g); %Number of triangles for each node
                %The local clustering coefficient of each node
                c = zeros(size(deg));
                c(deg > 1) = 2 * cn(deg > 1) ./ (deg(deg > 1).*(deg(deg > 1) - 1)); 
                %Average clustering coefficient of the graph
                %acc = mean(c(deg > 1)); 
                %acc = mean(c(deg > 1),'omitnan');
                ac(end+1) = mean(c);
                dummy = sum(cn) / sum(triu(g*g,1),'all');
                dummy(isnan(dummy)) = 0;
                gc(end+1) = dummy;
                
                steps(end+1) = j;
                
                disp(j)
            end
        end
    end
end