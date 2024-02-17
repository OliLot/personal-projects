disp(mean(tolsTiter))
disp(median(tolsTiter))

disp(mean(tolsBFGS))
disp(median(tolsBFGS))

%%%%%% ---------- Titer
figure
histogram(tolsTiter, 100, 'FaceAlpha', 0.5);

avg = mean(tolsTiter);
med = median(tolsTiter);
% xlim([0 1.5e-3]);
xline(avg,'--r', 'LineWidth', 1);
xline(med,'--m', 'LineWidth', 1);
xlabel('Tolerance Value');
ylabel('Bin Count');


figure
histogram(log10(tolsTiter), 50, 'FaceAlpha', 0.5);

avg = mean(log10(tolsTiter));
med = median(log10(tolsTiter));
%xlim([0 1.5e-3]);
xline(avg,'--r', 'LineWidth', 1);
xline(med,'--m', 'LineWidth', 1);
xlabel('Log Tolerance Value');
ylabel('Bin Count');


% -------------- BFGS
figure
histogram(tolsBFGS, 100, 'FaceAlpha', 0.5);

avg = mean(tolsBFGS);
med = median(tolsBFGS);
% xlim([0 1.5e-3]);
xline(avg,'--r', 'LineWidth', 1);
xline(med,'--m', 'LineWidth', 1);
xlabel('Tolerance Value');
ylabel('Bin Count');


figure
histogram(log10(tolsBFGS), 50, 'FaceAlpha', 0.5);

avg = mean(log10(tolsBFGS));
med = median(log10(tolsBFGS));
%xlim([0 1.5e-3]);
xline(avg,'--r', 'LineWidth', 1);
xline(med,'--m', 'LineWidth', 1);
xlabel('Log Tolerance Value');
ylabel('Bin Count');

% ---------- LARGE DF DISTRIBUTIONS for TITER

disp(mean(tolsLTiter))
disp(median(tolsLTiter))


figure
histogram(tolsLTiter, 100, 'FaceAlpha', 0.5);

avg = mean(tolsLTiter);
med = median(tolsTiter);
% xlim([0 1.5e-3]);
xline(avg,'--r', 'LineWidth', 1);
xline(med,'--m', 'LineWidth', 1);
xlabel('Tolerance Value');
ylabel('Bin Count');


figure
histogram(log10(tolsLTiter), 50, 'FaceAlpha', 0.5);

avg = mean(log10(tolsLTiter));
med = median(log10(tolsLTiter));
%xlim([0 1.5e-3]);
xline(avg,'--r', 'LineWidth', 1);
xline(med,'--m', 'LineWidth', 1);
xlabel('Log Tolerance Value');
ylabel('Bin Count');

