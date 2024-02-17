%{
Calculates the MLEs and their standard errors for all mixBvLap parameters 
given a Tx2 set of data.

INPUTS:
    - x = Tx2 dataset (T samples)
OUTPUT:
    - MLEs and their standard errors
%}
function [param, stderr, loglik] = mixBvLapMLE(x, initvec)
[n, d]=size(x); if d~=2, error('Data is not bivariate'), end

% Bounds for mu1_1, mu1_2, mu2_1, mu2_2, Sigma1_11, Sigma1_12&21, Sigma1_22,
% Sigma2_11, Sigma2_12&21, Sigma2_22, b1, b2, lam1
bound.lo= [-1 -1 -1 -1 0.01 -40 0.01 0.01 -50 0.01 0.01 0.01 0.01];
bound.hi= [1 1 1 1 40 40 40 40 40 40 40 40 1];
bound.which=[0 0 0 0 1 1 1 1 1 1 1 1 1];

% In this case, as bound.which for the mus are zero, mu will not be
% restricted. 
maxiter=500; tol=1e-6; MaxFunEvals=length(initvec)*maxiter;
opts=optimoptions('fminunc', 'Algorithm','quasi-newton','Display', ...
    'iter', 'HessianApproximation','bfgs','MaxIterations', ...
    maxiter,'OptimalityTolerance', tol, 'MaxFunctionEvaluations', MaxFunEvals);
[pout, fval, exitflag, output, grad, hess] = fminunc(@(param) ...
    mixBvLaploglik(param, x, bound), einschrk(initvec, bound), opts);

V=inv(hess)/n ; % Don't negate: we work with the neg of the loglik/n
[param, V]= einschrk(pout, bound, V) ; % Transform back, apply delta method
param=param';
loglik = -fval*n;
stderr=sqrt(diag(V)); % Approx std err of the params
end

%{
Calculates the negative log likelihood of a 2-component mixture of bivariate 
Laplace distributions given data x.

INPUTS:
    - x = Tx2 data set
    - param = array containing all distinct parameter values (16 total for
    4 location terms, 8 scale factors, 2 bs and 2 lams.)
OUTPUT:
    - log likelihood value
%}
function ll = mixBvLaploglik(param, x, bound)
if nargin < 3, bound=0; end
if isstruct(bound), paramvec=einschrk(real(param), bound, 999);
else paramvec=param;
end

% select parameters
mus1 = paramvec(1:2)';
mus2 = paramvec(3:4)';
Sigma1_11 = paramvec(5); Sigma1_12 = paramvec(6); Sigma1_22 = paramvec(7);
Sigma1_21 = Sigma1_12;
Sigma1 = [Sigma1_11 Sigma1_12; Sigma1_21 Sigma1_22];
Sigma2_11 = paramvec(8); Sigma2_12 = paramvec(9); Sigma2_22 = paramvec(10);
Sigma2_21 = Sigma2_12;
Sigma2 = [Sigma2_11 Sigma2_12; Sigma2_21 Sigma2_22];
b1 = paramvec(11); b2 = paramvec(12);
lam1 = paramvec(13); lam2 = 1-lam1;

% bivariate
d=2;

% Loop over density of each point
n = length(x(:,1));
if (min(eig(Sigma1)) < 1e-10) || (min(eig(Sigma2)) < 1e-10), ll=1e5; 
else
    pdfDensities = zeros(1, n);
    for i=1:n
        % disp(x(i,:));
        % disp(mus1);
        m1 = (x(i,:)' - mus1)' * Sigma1^(-1) * (x(i,:)' - mus1);
        m2 = (x(i,:)' - mus2)' * Sigma2^(-1) * (x(i,:)' - mus2);
        
        % get individual bivariate laplacian likelihoods
        pdfLap1 = (1 / (det(Sigma1)^(1/2) * (2*pi)^(d/2))) * 2/(gamma(b1)) * ...
            (m1/2)^(b1/2 - d/4) * besselk(b1 - d/2, sqrt(2*m1));
        pdfLap2 = (1 / (det(Sigma2)^(1/2) * (2*pi)^(d/2))) * 2/(gamma(b2)) * ... 
            (m2/2)^(b2/2 - d/4) * besselk(b2 - d/2, sqrt(2*m2));
        
        %mixture likelihood
        pdf_i = lam1*pdfLap1 + lam2*pdfLap2;
    
        %loglikelihood
        pdfDensities(i) = pdf_i;
    end
    logDensities = log(pdfDensities); ll = -mean(logDensities); 
    if isinf(ll), ll=1e5; end
end
end



% Einschrk
function [pout, Vout] = einschrk(pin, bound, Vin)
lo=bound.lo; hi=bound.hi; welche=bound.which;
if nargin < 3 %theta to phi
    trans=sqrt((hi-pin) ./ (pin-lo)); pout=(1-welche).* pin + welche .* trans;
    Vout =[];
else %phi to theta
    trans= (hi+lo .* pin.^2) ./ (1+ pin.^2); pout=(1-welche) .* pin + welche .* trans;
    % now adjust the standard errors
    trans=2*pin .* (lo-hi ) ./ (1+pin .^2) .^2; % partial derivatives of theta wrt phis
    d=(1-welche) + welche .* trans; % either unity or delta method .
    J=diag(d); Vout = J * Vin * J;
end
end