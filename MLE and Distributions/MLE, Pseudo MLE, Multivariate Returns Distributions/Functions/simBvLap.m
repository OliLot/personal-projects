%{
Function to generate n samples from a bivariate laplacian.

INPUTS:
    - n = number of samples
    - mus = 1x2 vector of the means
    - Sigma = scale matrix
    - b = gamma distribution scale factor
OUTPUT:
    - nx2 sample, where each row is a particular sample
%}
function x = simBvLap(n, mus, Sigma, b)
% input checks
% - dimensions
if ~isequal(size(mus), [1 2])
    error("Incorrect mus dimensions: should be 1x2"); 
end
if ~isequal(size(Sigma), [2 2]) 
    error("Incorrect Sigma dimensions: should be 2x2"); 
end
if length(b) ~= 1
    error("Incorrect bs dimensions: should be a scalar > 0"); 
end

% - values
if ~all(eig(Sigma) > 0)
    error("Sigma is not positive definite")
end
if b <= 0
    error("b parameter must be greater than zero")
end

% generate n bivariate samples
x = zeros(n, 2); % nx2 samples
g = gamrnd(b,1,1,n);
for i = 1:n
    samp = mvnrnd(mus, g(i)*Sigma);
    x(i,:) = samp; 
end
