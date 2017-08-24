%BACTERIOSENSITIVTYRUN Test the Peclet fitting and image analysis
%procedures by sampling a set of angles at a given Peclet number,
%generating simulated images for at the sampled angles, and then fitting a
%Peclet number to the resulting image fits.
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage

% Load marginal PDF and CDF
load('marginalPDF.mat')

% Sample nSamples angles at given Peclet number
Pec = 50;
nSamples = 180;
phiSample = SampleAngles(Pec,nSamples,marginalPDF);

% Generate simulated bacteriophage and fit
[simParams,fitParams] = SensitivityTest(phiSample);

% Calculate sample CDF
phi = linspace(-pi/2,pi/2,5e2);
sampCDF = discreteCDF(phiSample,phi);

% Fit Peclet number
fitPec = PecletFit(sampCDF,marginalCDF,phi);

% Save
save('simulatedResults.mat','simParams','fitParams','Pec','phiSample', ...
    'fitPec','-v7.3')


%%
function sampleCDF = discreteCDF(sample,phi)
sampleCDF = 0*phi';
for ii = 1:length(phi)
    dum = sample;
    dum(dum >= phi(ii)) = [];
    
    sampleCDF(ii) = length(dum);
end
% normalise
sampleCDF = sampleCDF / length(sample);
% spline
sampleCDF = csaps(phi,sampleCDF);
% sample
sampleCDF = fnval(sampleCDF,phi);

end