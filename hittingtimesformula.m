% Make two plots: one with the analytic hitting times as a function of s
% and one with the number of terms needed in the sum (for the analytic 
% hitting time) to reach machine epsilon as a function s.
% For each plot, set the target edit distance dtarget = 10 and use two 
% choices of the initial edit distance, d0 = 20 and d0 = 100.
% Save the first plot as 'hittingtimesformulah.fig' and
% 'hittingtimesformulah.eps'.
% Save the second plot as 'hittingtimesformulai.fig' and
% 'hittingtimesformulai.eps'.
% Dependencies: none

d0 = 20;
dtarget = 10;
s = 0:0.5:2000;
lens = length(s);
hs = zeros(1,lens);
is = zeros(1,lens);
for i = 1:lens
    [hs(i),is(i)] = hittingtime(d0,dtarget,s(i));
end

% PLOT
f1 = figure;
hold on
plot(s,hs)
%title('Expected hitting time')
xlabel('Rate parameter (s)')
ylabel('Expected hitting time')
set(gca,'fontsize',20)
%saveas(gcf,'hittingtimesformulah.eps','epsc')
%saveas(gcf,'hittingtimesformulah.fig')

f2 = figure;
hold on
plot(s,is)
%title('Expected hitting time')
xlabel('Rate parameter (s)')
ylabel('Number of terms needed')
set(gca,'fontsize',20)
%saveas(gcf,'hittingtimesformulai.eps','epsc')
%saveas(gcf,'hittingtimesformulai.fig')

d0 = 100;
dtarget = 10;
s = 0:0.5:2000;
lens = length(s);
hs = zeros(1,lens);
is = zeros(1,lens);
for i = 1:lens
    [hs(i),is(i)] = hittingtime(d0,dtarget,s(i));
end

figure(f1)
hold on
plot(s,hs)
%title('Expected hitting time')
xlabel('Rate parameter (s)')
ylabel('Expected hitting time')
set(gca,'fontsize',20)
legend('d_o = 20','d_o = 100','Location','northwest')
saveas(gcf,'hittingtimesformulah.eps','epsc')
saveas(gcf,'hittingtimesformulah.fig')

figure(f2)
hold on
plot(s,is)
%title('Expected hitting time')
xlabel('Rate parameter (s)')
ylabel('Number of terms needed')
set(gca,'fontsize',20)
legend('d_o = 20','d_o = 100','Location','northwest')
saveas(gcf,'hittingtimesformulai.eps','epsc')
saveas(gcf,'hittingtimesformulai.fig')

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