% Return adjacency matrix drawn from the SBM (stochastic block model).
% truth is a vector giving the desired labeling for each vertex.
% p is the in-cluster edge connection probability.
% q is the out-of-cluster edge connection probability.
function r = sbm(truth,p,q)
    n = length(truth);
    k = max(truth);
    r = zeros(n);

    indices = cell(1,k);
    lenindices = zeros(1,k);
    for i = 1:k
        indices{i} = find(truth == i);
        lenindices(i) = length(indices{i});
    end
    
    for i = 1:k
        r(indices{i}, indices{i}) = binornd(1,p,lenindices(i),lenindices(i));
        for j = i+1:k
            r(indices{i}, indices{j}) = binornd(1,q,lenindices(i),lenindices(j));
        end 
        
    end
    r = triu(r,1);
    r = r + r';
end