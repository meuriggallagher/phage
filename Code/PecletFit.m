%PECLETFIT Fits Peclet number by comparing sample and marginal CDFs
%   PE = PECLETFIT(SAMPLECDF,MARGCDF,PH) Fits the Peclet number by
%   comparing the SAMPLECDF and MARGCDF at angles PH through the use of
%   Kolmogorov-Smirnov
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function pe = PecletFit(sampleCDF,margCDF,ph)

global sampCDF marginalCDF phi
sampCDF = sampleCDF;
marginalCDF = margCDF;
phi = ph;

% NAG
nagParams = NagGlobSetup;

% Peclet limits
bl = 1e-2;
bu = 500;

[~, ~, ~, ~, ~, pe,~,~,~,~] = nag_glopt_bnd_mcs_solve( ...
    'MinFuncPecKS',nagParams{1,3}, bl, bu, nagParams{11,1},...
    nagParams{11,2},nagParams{11,3}, nagParams{1,1}, nagParams{1,2});

end