% 4.2 Computing MLE for 2-component mixture of bivariate laplace
% distribution via BFGS, also reporting approximate standard errors.
rng(1, 'twister')

n = 10000; % sample size

x = simMixBvLap(n, [0 0], [0 0], [1.5 0.5; 0.5 1], [20 0; 0 10], [10 5], [0.9 0.1]);

initvec = [0 0 0 0 1 0.1 1 10 0 10 10 5 0.8]; % reasonable values extractable from the data; 10 5 0.8 are standard.
[param, stderr, loglik] = mixBvLapMLE(x, initvec); %1:2=mus1, 3:4=mus2, 5:7=Sig1, 8:10=Sig2, 11:12 = bs, 13=lam1
