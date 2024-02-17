%{
Function to simulate from a 2-component mixed bivariate laplacian
distribution.

INPUTS:
    - n = sample size
    - mus1 = 1x2 vector of means for the first component
    - mus2 = 1x2 vector of means for the second component
    - Sigma1 = 2x2 scale matrix for the first component
    - Sigma2 = 2x2 scale matrix for the second component
    - bs = 1x2 vector of gamma factors for the first and second component
    - lam = 1x2 vector of the mixture weights
OUTPUT:
    - nx2 sample from the mixed laplacian distribution
%}
function y = simMixBvLap(n, mus1, mus2, Sigma1, Sigma2, bs, lam)
% check inputs
% - dimensions
if ~isequal(size(mus1), [1 2]) % mean vector
    error("Incorrect mus1 dimensions: should be 1x2"); 
end
if ~isequal(size(mus2), [1 2]) % mean vector
    error("Incorrect mus2 dimensions: should be 1x2"); 
end
if ~isequal(size(Sigma1), [2 2]) % scale matrix
    error("Incorrect Sigma1 dimensions: should be 2x2"); 
end
if ~isequal(size(Sigma2), [2 2]) % scale matrix
    error("Incorrect Sigma2 dimensions: should be 2x2"); 
end
if ~isequal(size(bs), [1 2]) %
    error("Incorrect bs dimensions: should be 1x2"); 
end
if ~isequal(size(lam), [1 2]) % mixture weights
    error("Incorrect mixture weights dimension: should be 1x2"); 
end
% - values
if (any(lam < 0)) || (any(lam > 1)) || (abs(sum(lam) - 1) >= 1e-10) 
    error("Incorrect weight inputs: Either values less than zero, greater than one" + ...
        " or do not sum up to one.")
end
if ~all(eig(Sigma1) > 0), error("Sigma1 not positive definite"); end
if ~all(eig(Sigma2) > 0), error("Sigma2 not positive definite"); end
if any(bs <= 0), error("Incorrect b parameters: must be greater than zero."); end


k = 2; % k=2 mixture components
pick = zeros(n, k); % n samples, k mixtures 

% simulate realisations of multinomial for picking which mixture to use
for i = 1:n
    mult = randmultinomial(lam); 
    pick(i, mult) = 1;
end

% split into separate columns and double up as the samples are bivariate
pick1 = pick(:, 1);
pick1 = [pick1 pick1]; 
pick2 = pick(:, 2);
pick2 = [pick2 pick2];

% get both mixture samples
L1 = simBvLap(n, mus1, Sigma1, bs(1));
L2 = simBvLap(n, mus2, Sigma2, bs(2));

% implement mixing
y = pick1.* L1 + pick2.* L2; % select & combine realisations
end

% multinomial selection function
function C = randmultinomial(p), pp=cumsum(p); u=rand; C=1+sum(u>pp); end