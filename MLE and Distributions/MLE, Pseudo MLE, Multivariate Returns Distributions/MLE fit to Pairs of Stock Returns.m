% Task 4.4: randomly sample two stocks and fit both ML and NCT, collect AIC
% and BIC values and plot results.

% Get data
data = load("DJIA30stockreturns.mat").DJIARet;

% run function
tableAICBIC = fitMLandNCT(data, 50);
disp(tableAICBIC)
