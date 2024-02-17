function [timeTiter, timeBFGS, tolsTiter, tolsBFGS] = mleTimes(df, loc, sig, T, iter)

if nargin == 1
    T = 250;
    iter = 1000;
    loc = 0;
    sig = 1;
end
if nargin == 3
    T = 250;
    iter = 1000;
end

% 1. timing titer
rng(1, "twister")
timesTiter = zeros(1, iter);
tolsTiter = zeros(1, iter);
for i=1:iter
    % generate dataset
    x = loc + sig*trnd(df, T, 1); % generalised t distribution sample
    
    % find 4 sf tolerances
    tolTiter = tolFinder(@titer, x);
    tolsTiter(i) = tolTiter;
    
    % time calculations
    tic
    [dfhat, muhat, chat, iters1] = titer(x, tolTiter, 1);
    timesTiter(i) = toc;
end
timeTiter = sum(timesTiter);

% 2. timing BFGS
rng(1, "twister")
timesBFGS = zeros(1, iter);
tolsBFGS = zeros(1, iter);
for i=1:iter
    % generate dataset
    x = loc + sig*trnd(df, T, 1); % generalised t distribution sample
    
    % find 4 sf tolerances
    tolBFGS = tolFinder(@tlikmax, x, [2, 2, 2]);
    tolsBFGS(i) = tolBFGS;
    
    % time MLE
    tic
    MLE = tlikmax(x, tolBFGS, [2, 2, 2]);
    timesBFGS(i) = toc;
end
timeBFGS = sum(timesBFGS);