% Compute the analytic expected hitting time given d0 (initial edit
% distance), dtarget (target edit distance), and s (rate parameter).
% Return h (the expected hitting time) and i (the index of the first term 
% in the sum smaller than machine epsilon).
function [h,i] = hittingtime(d0,dtarget,s)
    term = exp(-1/s) * (1 - exp(-(d0-dtarget)/s))/(1 - exp(-1/s));
    sums = term;
    i = 1;
    % Assume dm (the maximum edit distance) is sufficiently large so the 
    % number of terms in the sum is limited by machine precision.
    while term > eps
        i = i + 1;
        term = exp(-sum(1:i)/s) * (1 - exp(-i*(d0-dtarget)/s)) / (1 - exp(-i/s));
        sums = sums + term;
    end
    h = d0 - dtarget + 2 * sums;
end