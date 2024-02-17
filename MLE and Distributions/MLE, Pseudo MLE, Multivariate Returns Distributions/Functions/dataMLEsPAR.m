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
function [paramML, aicML, bicML, paramNCT, aicNCT, bicNCT] = dataMLEsPAR(x)
% pos mean, pos corr
% pos mean, neg corr
% neg mean, pos corr
% neg mean, neg corr
% pos/neg mean, pos corr
%pos/neg mean, neg corr
% neg/pos mean, pos corr
% neg/pos mean, neg corr

[n, d] = size(x); if d ~= 2, error("data is not nx2"); end

allParamML = zeros(13, 2);
allStdErr = zeros(13, 2);
AllLogLik = zeros(2);
AllInitialVectors = [
    [0.1 0.1 0.1 0.1 2 1 4 10 1 10 10 5 0.9];
    [0.1 0.1 0.1 0.1 4 -1 2 10 -1 10 10 5 0.9];
    ];

parfor i = 1:2
initVec = AllInitialVectors(i, :)
[paramML, stderrML, loglik] = mixBvLapMLE(x, initVec);
allParamML(:, i) = paramML; allStdErr(:, i) = stderrML; AllLogLik(i) = loglik;
end

%params = [paramML_pp paramML_pn paramML_np paramML_nn paramML_pnp paramML_pnn paramML_npp paramML_npn];
%errs = [stderrML_pp stderrML_pn stderrML_np stderrML_nn stderrML_pnp stderrML_pnn stderrML_npp stderrML_npn];
%loglikelihoodsML = [loglik_pp, loglik_pn, loglik_np, loglik_nn, loglik_pnp, loglik_pnn, loglik_npp, loglik_npn];

params = allParamML;
errs = allStdErr;
loglikelihoodsML = AllLogLik;

for i = 1:2 % 2 = number of tested initialisations
    if isreal(errs(:, i)), continue
    elseif isnan(errs(:, i)), errs(:, i) = 1000*ones(13, 1);
    else errs(:, i) = 1000*ones(13, 1); % 13 = number of parameters
    end 
end

sqSumErrs = sum(errs.^2);
[~, indParML] = min(sqSumErrs);
disp(indParML)

paramML = params(:,indParML);
loglikML = loglikelihoodsML(indParML);


[paramNCT, stderrNCT, loglikNCT] = bvNCTMLE(x);

[aicML, bicML] = aicbic(loglikML, length(paramML), n);
[aicNCT, bicNCT] = aicbic(loglikNCT, length(paramNCT), n);
