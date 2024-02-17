%get table values
aicsML = tableAICBIC.aicsML;
aicsNCT = tableAICBIC.aicsNCT;
bicsML = tableAICBIC.bicsML;
bicsNCT = tableAICBIC.bicsNCT;

% proportions preferred
aicPrefML = sum(aicsML < aicsNCT) / length(aicsML);
disp(aicPrefML)
bicPrefML = sum(bicsML < bicsNCT) / length(bicsML);
disp(bicPrefML)

% plots
figure
scatter(aicsML', aicsNCT')
xlim([4.5e4 7e4])
ylim([4.5e4 7e4])
hold on
plot(xlim,ylim,'-k')
xlabel('aicML')
ylabel('aicNCT')
title('Model AICs')
hold off

figure
scatter(bicsML', bicsNCT')
xlim([4.5e4 7e4])
ylim([4.5e4 7e4])
hold on
plot(xlim,ylim,'-k')
xlabel('bicML')
ylabel('bicNCT')
title('Model BICs')
hold off