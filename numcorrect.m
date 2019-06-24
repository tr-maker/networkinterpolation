% Given the labeling X and the true labeling truthinput, return the number
% of correct labels in X and the permutation of the true labeling that best
% aligns with X.

function [corr,truth] = numcorrect(X,truthinput)
    n = numel(X);
    if n ~= numel(truthinput)
        error(['The given labeling and true labeling have different '...
            'numbers of elements']);
    end
    
    truth = truthinput;
    corr = length(find(X == truth));
    if corr == n
        return;
    end

    k = max(truthinput);
    idx = cell(1,k); %index set of each label in the truthinput
    for i = 1:k
        idx{i} = find(truthinput == i);
    end
    a = nextperm(1:k);

    while (~isempty(a))
        newtruth = zeros(1,n);
        for i = 1:k
            newtruth(idx{i}) = a(i);
        end

        newcorr = length(find(X == newtruth));
        if newcorr > corr
            corr = newcorr;
            truth = newtruth;
            if corr == n
                return;
            end
        end
        a = nextperm(a);
    end
    return;
end

% Return the next permutation (in lexicographic order) of the given row
% permutation vector
function b = nextperm(a)
    k = length(a);
    i = find(a(2:k) > a(1:k-1), 1, 'last');
    if isempty(i)
        b = [];
    else
        j = find(a > a(i), 1, 'last');
        a([i j]) = a([j i]);
        a(i+1:k) = fliplr(a(i+1:k));
        b = a;
    end 
end