%{
Function computing MLE and Standard errors for the bivariate NCT
distribution.

INPUT:
    - x = Tx2 dataset, where T is the number of samples
OUTPUT:
    - param = MLEs
    - stderr = approximation of the standard errors
    - iters = number of iterations
    - loglik = loglikelihood at final value
%}
function [param,stderr,loglik] = bvNCTMLE(x)
x = x';
[d, T]=size(x); if d~=2, error('Data is not bivariate'), end
%%%%%%%% k R12 gam1 gam2
bound.lo= [ 1.1 -1 -4 -4 ];
bound.hi= [ 20 1 4 4 ];
bound.which=[ 1 1 1 1 ];
initvec =[ 2 0.5 1 1 ];


maxiter=300; tol=1e-6; MaxFunEvals=length(initvec)*maxiter;
opts=optimset('Display','iter','Maxiter',maxiter,'TolFun',tol,'TolX',tol, ...
    'MaxFunEvals',MaxFunEvals,'LargeScale','Off');
[pout,fval,~,~,~,hess]= ...
    fminunc(@(param) MVNCTloglik(param,x,bound),einschrk(initvec,bound),opts);
V=inv(hess)/T; [param,V]=einschrk(pout,bound,V); param=param'; 
stderr=sqrt(diag(V)); loglik=-fval*T;
end


function ll=MVNCTloglik(param,x,bound)
if nargin<3, bound=0; end
if isstruct(bound), param=einschrk(real(param),bound,999); end
k=param(1); R12=param(2); gam=param(3:4); 
R=[1 R12; R12 1]; 
if min(eig(R))<1e-4, ll=1e5; 
else
    llvec = mvnctpdfln(x, gam, k, R); 
    ll=-mean(llvec); if isinf(ll), ll=1e5; end
end

end

function pdfln = mvnctpdfln(x, gam, v, Sigma)
% x = Tx2 matrix of evaluation points
% gam = 2-length noncentrality vector
% v = degrees of freedom
% Sigma = dispersion matrix

% input checks
if (~isequal(size(v), [1 1])) || (v <= 0)
    error("Incorrect v dimensions: should be a scalar greater than one"); 
end
if ~isequal(size(gam), [1 2])
    error("Incorrect gam dimensions: should be 1x2"); 
end
if ~isequal(size(Sigma), [2 2])
    error("Incorrect Sigma dimensions: should be 2x2 "); 
end
if ~all(eig(Sigma) >= 0), error("Sigma not positive semi-definite"); end

% main body of function
[d,~] = size(x); C=Sigma; [R, ~] = cholcov(C, 0);
gam=reshape(gam,length(gam),1);
vn2 = (v + d) / 2; xm = x; rho = sum((R'\xm).^2,1);
pdfln = gammaln(vn2) - d/2*log(pi*v) - gammaln(v/2) - ...
sum(slog(diag(R))) - vn2*log1p(rho/v);
if (all(gam == 0)), return; end
idx = (pdfln >= -37); maxiter=1e4; k=0;
if (any(idx))
    gcg = sum((R'\gam).^2); pdfln = pdfln - 0.5*gcg; xcg = xm' * (C \ gam);
    term = 0.5*log(2) + log(xcg) - 0.5*slog(v+rho');
    term(term == -inf) = log(realmin); term(term == +inf) = log(realmax);
    logterms = gammaln((v+d+k)/2) - gammaln(k+1) - gammaln(vn2) + k*term;
    ff = real(exp(logterms)); logsumk = log(ff);
    while (k < maxiter)
        k=k+1;
        logterms = gammaln((v+d+k)/2) - gammaln(k+1) - gammaln(vn2) + k*term(idx);
        ff = real(exp(logterms-logsumk(idx))); logsumk(idx)=logsumk(idx)+log1p(ff);
        idx(idx) = (abs(ff) > 1e-4); if (all(idx == false)), break, end
    end
    pdfln = real(pdfln+logsumk');
end
end


function y = slog(x) % Truncated log. No -Inf or +Inf.
y = log(max(realmin, min(realmax, x)));
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
