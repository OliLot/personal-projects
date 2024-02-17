function sample = simBvNCT(n, v, y, R)
sample = zeros(n,2);

for i=1:n
    Z = mvnrnd(y, R);
    C = chi2rnd(v, 1);
    sample(i, :) = Z / (sqrt(C/v));
end
