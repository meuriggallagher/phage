%SAMPLEANGLES Samples angles from the marginal PDF
%   PHISAMPLES = SAMPLEANGLES(PEC,NSAMPLES,MARGINALPDF) Samples NSAMPLES
%   angles at Peclet number PEC from MARGINALPDF
%
%   INPUTS:
%       PEC         - Peclet number to sample at
%       NSAMPLES    - Number of samples to take
%       MARGINALPDF - Marginal PDF from which to sample
%
%   OUTPUTS:
%       PHISAMPLES  - Sampled angles
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function phiSamples = SampleAngles(Pec,nSamples,marginalPDF)

%% Calculate CDF for given peclet number
nS = 100;

phi = linspace(-pi/2,pi/2,2*nS-1);
dP = phi(3) - phi(1);

% Sample PDF
PDF = fnval(marginalPDF,[phi;Pec*ones(1,length(phi))]);
PDF(PDF < 0) = 1e-20;

% Integrate to find CDF
CDF = zeros(nS,1);
for ii = 1:nS-1
    CDF(ii+1) = dP/6 * (PDF(2*ii-1) + 4*PDF(2*ii) + PDF(2*ii+1));
end
CDF = cumsum(CDF);

% Scale to 1 at end
CDF = CDF * 1/CDF(end);

% Spline CDF
phi = linspace(-pi/2,pi/2,nS);
CDFspline = csapi(phi,CDF);

%% Sample points from CDF
% Use INVERSE TRANSFORM SAMPLING to generate pseudo-random samples
samples = rand(nSamples,1);
phiSamples = zeros(nSamples,1);
for ii = 1:length(samples)

    CDFfunc = @(p) fnval(CDFspline,p) - samples(ii);
    
    phiSamples(ii) = fzero(CDFfunc,0);
end

end