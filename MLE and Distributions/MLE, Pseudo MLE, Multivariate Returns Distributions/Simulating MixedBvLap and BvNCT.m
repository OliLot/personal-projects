% 4.1 Simulate 2-component mixture of bivariate laplace distributions
rng(1,'twister')

% stocks with positive mean return and correlation 1/9, extreme events have
% mean zero and slight correlation. Extreme events have a probability of 
% 0.05.
x = simMixBvLap(100000, [0 0], [0 0], [1.5 0.5; 0.5 1], ...
    [20 0; 0 10], [10 5], [0.9 0.1]);

figure;
hist3(x, 'ctrs',{-20:0.5:20 -20:0.5:20});

y = simBvNCT(100000, 3, [0.4 -0.3], [1 0.3; 0.3 1]);

figure;
hist3(y, 'ctrs',{-10:0.25:10 -10:0.25:10});

