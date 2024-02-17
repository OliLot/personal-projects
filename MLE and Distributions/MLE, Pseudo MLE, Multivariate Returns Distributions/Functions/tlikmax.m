function [param, stderr, iters, loglik, Varcov] = tlikmax(x, tol, initvec)

%%%%%%%% Getting df, mu, c
bound.lo= [1 1 0.01];
bound.hi= [100 1 100];
bound.which=[1 0 1];

% In this case, as bound.which for mu is zero, mu will not be
% restricted. As such, the values for .lo and  hi are irrelevant
maxiter = 200; % change these as you see fit
opts=optimoptions('fminunc', 'Algorithm','quasi-newton','Display','notify-detailed', ...
    'HessianApproximation','bfgs','MaxIterations', maxiter,'OptimalityTolerance', tol);
[pout, fval, exitflag, theoutput, grad, hess] = fminunc(@(param) tloglik(param, x, bound), einschrk(initvec, bound), opts);

V=inv(hess) ; % Don't negate: we work with the neg of the loglik
[param, V]= einschrk(pout, bound, V) ; % Transform back , apply delta method
param=param'; Varcov=V;
stderr=sqrt(diag(V)); % Approx std err of the params
loglik=-fval; % The value of the loglik at its maximum.
iters=theoutput.iterations; % Number of loglikfunction evals


% Loglikelihood
function ll = tloglik (param, x, bound)
if nargin < 3, bound=0; end
if isstruct(bound), paramvec=einschrk(real(param), bound, 999);
else paramvec=param;
end

v=paramvec(1); mu=paramvec(2); c=paramvec(3);
K=beta(v/2, 0.5) * sqrt(v); z=(x-mu)/c;
ll= -sum(-log(c) - log(K) - ((v+1)/2)*log(1 + (z.^2)/v));



% Einschrk
function [pout, Vout] = einschrk(pin, bound, Vin)
lo=bound.lo; hi=bound.hi; welche=bound.which;
if nargin < 3
    trans=sqrt((hi-pin) ./ (pin-lo)); pout=(1-welche).* pin + welche .* trans;
    Vout =[];
else
    trans= (hi+lo .* pin.^2) ./ (1+ pin.^2); pout=(1-welche) .* pin + welche .* trans;
    % now adjust the standard errors
    trans=2*pin .* (lo-hi ) ./ (1+pin .^2) .^2;
    d=(1-welche) + welche .* trans; % either unity or delta method .
    J=diag(d); Vout = J * Vin * J;
end
