%{
Function running MLE fits on two random data samples from x for the
2-component bivariate Laplacian mixture and the NCT, outputting the AIC and
BIC values of each fit. This is repeated rep number of times, with the
function outputting arrays of the AIC and BIC values (column 1 is ML and
comlumn 2 is NCT) as well as plotting the AIC and BIC results.

INPUTS:
    - x = N x b data matrix where N is the number of samples and b are the
    different data vectors
    - rep = number of MLE fit repetitions to run
OUTPUTS:
    - AICs = rep x 2 array for ML and NCT MLE fit AIC values
    - BICs = rep x 2 array for ML and NCT MLE fit BIC values
    - two plots showing the AIC and BIC trends
%}
function [tableAICBICs] = fitMLandNCT(x, rep)

% get the total number of data vectors
[n, b] = size(x);

% initialise AIC and BIC arrays
AICBICs = zeros(rep, 4);
% repeat random 2 samples from data
for i=1:rep
    % sample indices
    [ind1, ind2] = randPair(1, b);
    
    disp([ind1 ind2])

    % sample data
    x1 = x(:, ind1);
    x2 = x(:, ind2);
    dat = [x1 x2];

    % BvML MLE and the NCT MLE
    [~, aicML, bicML, ~, aicNCT, bicNCT] = dataMLEsPAR(dat); % do not need parameters

    % input
    AICBICs(i, 1:2) = [aicML aicNCT];
    AICBICs(i, 3:4) = [bicML bicNCT];
end
aicsML = AICBICs(:, 1);
aicsNCT = AICBICs(:, 2);
bicsML = AICBICs(:, 3);
bicsNCT = AICBICs(:, 4);

% output table
tableAICBICs = table(aicsML, aicsNCT, bicsML, bicsNCT);



% Functions plotting the AIC and BIC values
figure
scatter(aicsML', aicsNCT')
plot(xlim,ylim,'-b')
xlabel('aicML')
ylabel('aicNCT')
title('Model AICs')

figure
scatter(bicsML', bicsNCT')
plot(xlim,ylim,'-b')
xlabel('bicML')
ylabel('bicNCT')
title('Model BICs')

end



% random pair of integers between a and b
function [ind1, ind2] = randPair(a, b)
ind1 = randi([a b], 1);
ind2 = randi([a b], 1);
while ind2 == ind1
    ind2 = randi([a b], 1);
end
end