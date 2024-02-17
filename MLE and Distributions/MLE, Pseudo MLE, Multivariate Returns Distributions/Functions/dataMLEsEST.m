%{
Function that, for a given dataset, computes the MLEs of a 2-component
mixture of bivariate Laplacian distributions and of a bivariate NCT
distribution.

INPUT:
    - x = Tx2 dataset where T is the sample size
OUTPUT:
    - paramML = MLEs for the mixed laplace distribution
    - stderrML = approximate standard errors for the mixed laplace MLEs
    - paramNCT = MLEs for the bivariate NCT distribution
    - stderrML = approximate standard errors for the bivariate NCT MLEs
%}
function [paramML, aicML, bicML, paramNCT, aicNCT, bicNCT] = dataMLEsEST(x)

[n, d] = size(x); if d~=2, error("incorrect number of columns, should be 2"); end

% Get rough estimates for initial vectors
x1 = x(:, 1)';
mean1 = mean(x1);
med1 = median(x1);
std1 = std(x1);
x1Sort = sort(x1);
std1Sort = std(x1Sort(n/10:9*n/10)); %cut out 1/10 th of the lowest and highest values -WORKS WELL

x2 = x(:, 2)';
mean2 = mean(x2);
med2 = median(x2);
std2 = std(x2);
x2Sort = sort(x2);
std2Sort = std(x2Sort(n/10:9*n/10)); %cut out 1/10 th of the lowest and highest values -WORKS WELL

Cov12 = cov(x1, x2);
covVal = Cov12(1, 2);
covVal = sign(covVal); % set equal to 1 or -1 depending on predicted relationship


initvecML = [med1 med2 mean1 mean2 std1Sort/2 covVal std2Sort/2 2*std1 covVal 2*std2 10 5 0.8]; % 10, 5, 0.8 standard values
disp(initvecML)
[paramML, ~, loglikML] = mixBvLapMLE(x, initvecML); 

% x = simMixBvLap(10000, paramML(1:2), paramML(3:4), [paramML(5) paramML(6); ... 
%     paramML(6) paramML(7)], [paramML(8) paramML(9); paramML(9) paramML(10)], ...
%     paramML(11:12), [paramML(13) 1 - paramML(13)]); % n mu1 mu2 S1 S2 bs lams
% 
% figure
% hist3(x, 'nbins', [100 100])

% only those with real stderrs; print these; from these take smallest sum
% of squared errors as the answer?

[paramNCT, ~, loglikNCT] = bvNCTMLE(x); % k R12 gam1 gam2

% y = simBvNCT(10000, paramNCT(1), paramNCT(3:4), [1 paramNCT(2); paramNCT(2) 1]); % n k gam R
% hist3(y, 'nbins', [100 100])

% get AIC and BIC values
[aicML, bicML] = aicbic(loglikML, length(paramML), n);
[aicNCT, bicNCT] = aicbic(loglikNCT, length(paramNCT), n);
