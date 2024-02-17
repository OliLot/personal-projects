%4.3 BivariateNCT MLE results
rng(1, 'twister')

x = simMixBvLap(1000, [10 20], [-10 -20], [3 2; 2 5], [13 2; 2 20], [10 5], [0.9 0.1]);
hist3(x,'nbins',[100 100]);

y = simBvNCT(1000, 4, [2 2], [1 0.5; 0.5 1]);
hist3(y, 'nbins', [100 100])

% [xparamML, xindParML, xloglikML, xparamNCT, xloglikNCT] = dataMLEsPAR(x);
% [yparamML, yindParML, yloglikML, yparamNCT, yloglikNCT] = dataMLEsPAR(y);

[xparamML, xaicML, xbicML, xparamNCT, xaicNCT, xbicNCT] = dataMLEsEST(x);
[yparamML, yaicML, ybicML, yparamNCT, yaicNCT, ybicNCT] = dataMLEsEST(y);


% k mu1 mu2 scale1 scale2 R12 gam1 gam2
disp(xparamML)
disp(xparamNCT)
disp([xaicML, xbicML, xaicNCT, xbicNCT])

disp(yparamML)
disp(yparamNCT)
disp([yaicML, ybicML, yaicNCT, ybicNCT])