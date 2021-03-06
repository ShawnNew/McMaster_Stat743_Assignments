clear;clc;
%% Initialization

n = 1:1:100;   % test the shorter CI is not change as n varies
alpha = 0.05; % for a 95% confidence CI

%% Define two functions of the two lengths
L_pivot = @(alpha, n) (alpha/2)^(-1/n) - (1-alpha/2)^(-1/n);
L_LRT = @(alpha, n) alpha^(-1/n) - 1;

len = zeros(1, length(n));

for i = 1 : length(n)
    %len(1, i) = L_pivot(alpha, i);
    %len(2, i) = L_LRT(alpha, i);
    len(i) = L_pivot(alpha, i) / L_LRT(alpha, i);
end

%% Infinity
syms i
lim = limit(L_pivot(alpha, i)/L_LRT(alpha, i), i, Inf);


%% Draw the picture
plot(n, len, 'b--o')
grid on;
title('The ratio of two lengths(pivotal/LRT) of the interval');
xlabel('n');
ylabel('ratio');

