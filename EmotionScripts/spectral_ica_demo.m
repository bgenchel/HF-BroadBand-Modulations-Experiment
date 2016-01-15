%% calculate spectral decomposition for a subset of components

comps = [17 25 5 6 7 8 9 11 12 13 15 16 17 18 19 21 22 27];

clear icSpectralDecomposition imActivityEpochs logPowerModulationEpochs powerForWavelet powerForWaveletEpochs 
for compId = 1:length(comps)
    compId
    icSpectralDecomposition(compId) = make_spectral_decomposition(ALLEEG(1), comps(compId), 10);
end;

%% plotting the template

% plot all the templates
plot_spectral_template(icSpectralDecomposition(1)); 

% plot a single template 
figure;
plot_single_spectral_template(icSpectralDecomposition(1),5);
%% calculate spectal templte activity correlation and show in a figure

combinedTemplateActivation = [];
for i=1:length(icSpectralDecomposition)
    cleanActivityI = icSpectralDecomposition(i).activity;
    cleanActivityI(:,icSpectralDecomposition(i).badFrames) = [];
    combinedTemplateActivation = cat(1, combinedTemplateActivation, cleanActivityI);
end;

combinedTemplateActivationCorrelation = corrcoef(combinedTemplateActivation');

% remove ones from the diagonal of the correlation matrix (for better visualization).
combinedTemplateActivationCorrelation = combinedTemplateActivationCorrelation -diag(diag(combinedTemplateActivationCorrelation));


figure;
imagesc((combinedTemplateActivationCorrelation));
title('Correlation between activity of spectral modes of cortical ICs','fontsize',16);
