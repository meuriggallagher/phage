%SENSITIVITYTEST Simulated phage fit
%   [SIMPARAMS,FITPARAMS] = SENSITIVITYTEST(PHI) Generates a simulated
%   bacteriophage image at angles PHI, with parameters SIMPARAMS and fits
%   to the resulting image giving parameters FITPARAMS
%
%   INPUTS:
%       PHI       - Set of bacteriophage angles to simulate
%
%   OUTPUTS:
%       SIMPARAMS - Simulated parameters for bacteriophage
%       [x0;y0;phi;rho;sX;background;error]
%       FITPARAMS - Fit parameters from simulated image
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function [simParams,fitParams] = SensitivityTest(phi)

% Memory allocation
nPhi = length(phi);
simParams = zeros(7,nPhi);
fitParams = zeros(7,nPhi);

% Image size
nX = 110; nY = 110;
[~,~,~,x,y] = ImageParameters(nX,nY);

% Fixed variables
rho = 1.5e-6;
sX = ImagingParameters*3;
nS = 35;
s = linspace(0,rho,2*nS-1);
backgroundStd = 5;

% Simulated parameters
simParams(3,:) = phi;
simParams(4,:) = rho*ones(nPhi,1);
simParams(5,:) = sX*ones(nPhi,1);
simParams(7,:) = backgroundStd*ones(nPhi,1);
backgroundFactor = 0.3;

% Generate image for each bacteriophage
for ii = 1:size(simParams,2)
    fprintf('\tFrame %i of %i\n',ii,nPhi)
    
    % Random parameters
    x0Lim = [-1,1]*1e-6;
    y0Lim = x0Lim;
    x0 = x0Lim(1) + (x0Lim(2)-x0Lim(1)).*rand;
    y0 = y0Lim(1) + (y0Lim(2)-y0Lim(1)).*rand;
    
    % Generate image
    simImage = GenImage(x,y,x0,y0,phi(ii),s,sX,nS,0);
    
    % Scale simimulated image to [0,255]
    simImage = simImage.*(255/max(simImage(:)));
    
    % background
    backgroundMean = max(simImage(:)) * backgroundFactor;
    simBkgrnd = backgroundMean + backgroundStd*randn(nY,nX);
    simImage = simImage + simBkgrnd;
    
    % round image
    simImage = round(simImage); simImage(simImage < 0) = 0;
    
    % Preprocess image
    imageFrame = simImage;
    imageFrame = medfilt2(imageFrame,[5 5]);
    imageFrame = imageFrame - median(imageFrame(:));
    imageFrame(imageFrame < 3) = 0;
    
    % Sim params
    simParams(1,ii) = x0;
    simParams(2,ii) = y0;
    simParams(6,ii) = backgroundMean;
    
    tic
    fitParams(:,ii) = Phage2DFit(imageFrame);
    t1 = toc;
    fprintf('\t\tTime taken = %f\n',t1)
    fprintf('\tsimulated phi = %e, fit phi = %e\n',simParams(3,ii), ...
        mod(fitParams(3,ii)+pi/2,pi)-pi/2)
    fprintf('\n\tDiff = %e\n\n',simParams(3,ii) - ...
        (mod(fitParams(3,ii)+pi/2,pi)-pi/2))
    
end