clear;clc;
% Initialization
t = 0.85;
n = 1:1:100;   % test the shorter CI is not change as n varies
alpha = 0.05; % for a 95% confidence CI

% Define two functions of the two lengths
L_pivot = @(t, alpha, n) t * (1/((alpha/2)^(1/n)) - 1/((1-alpha/2)^(1/n)));
L_LRT = @(t, alpha, n) t * (1/(alpha^(1/n)) - 1);

len = zeros(2, length(n));

for i = 1 : length(n)
    len(1, i) = L_pivot(t, alpha, i);
    len(2, i) = L_LRT(t, alpha, i);
end

% Draw the picture
plot(n, log(len(1,:)), 'b--o', n, log(len(2,:)), 'c--*')
grid on;
title('The difference between two lengths of the interval (suppose t is 0.85)');
xlabel('n');
ylabel('length');
legend('Length for pivot method', 'Length for LRT method');
