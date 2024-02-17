% // Part 1: Time comparison for the same tolerance of titer and MLE
% generate 1 year (250 days) of random student t data


% testing location effect on time: df = 5, mu = 0, sig = 1
[timeTiter, timeBFGS, tolsTiter, tolsBFGS] = mleTimes(4);
times = [timeTiter timeBFGS]';

% testing location effect on time: df = 15, mu = 0, sig = 1
[timeLTiter, timeLBFGS, tolsLTiter, tolsLBFGS] = mleTimes(15);
timesL = [timeLTiter timeLBFGS]';

