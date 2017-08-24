%MINFUNCOPTS Objective function for fitting the PSF width
%   [F,USER,INFORM]=MINFUNCOPTS(N,PARAMETERS,NSTATE,USER Objective function
%   for fitting the Gaussian PSF width parameter, set up for the NAG global 
%   optimisation regime e05jbc
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function [f,user,inform] = MinFuncOpts(n,parameters,nstate,user) ...
    %#ok<INUSL>

global imageFrame x y x0 y0 rho phi background nS

sX_1 = parameters(1);

s = linspace(0,rho,2*nS - 1);

obj = GenImage(x,y,x0,y0,phi,s,sX_1,nS,background);
obj = obj - background;
obj = obj ./ max(obj(:)) * (max(imageFrame(:))-background);
obj = obj + background;
f = sum((imageFrame(:)-obj(:)).^2);


inform = int64(1);

end