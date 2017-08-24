%MINFUNCINTS Objective function for fitting the background 
%   [F,USER,INFORM]=MINFUNCINTS(N,PARAMETERS,NSTATE,USER Objective function
%   for fitting the background parameter, set up for the NAG global
%   optimisation regime e05jbc
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function [f,user,inform] = MinFuncInts(n,parameters,nstate,user) ...
    %#ok<INUSL>

global imageFrame x y x0 y0 rho phi sX nS

background_1 = parameters(1);

s = linspace(0,rho,2*nS - 1);

obj = GenImage(x,y,x0,y0,phi,s,sX,nS,background_1);
obj = obj - background_1;
obj = obj ./ max(obj(:)) * (max(imageFrame(:))-background_1);
obj = obj + background_1;
f = sum((imageFrame(:)-obj(:)).^2);

inform = int64(1);

end