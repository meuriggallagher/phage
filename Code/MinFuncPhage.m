%MINFUNCPHAGE Objective function for fitting the bacteriophage parameters
%   [F,USER,INFORM]=MINFUNCPHAGE(N,PARAMETERS,NSTATE,USER Objective function
%   for fitting the bacteriophage parameters (x0,y0,rho,phi), set up for 
%   the NAG global optimisation regime e05jbc
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function [f,user,inform] = MinFuncPhage(n,parameters,nstate,user) ...
    %#ok<INUSL>

global imageFrame x y sX background nS

x0_1 = parameters(1);
y0_1 = parameters(2);
rho_1 = parameters(3);
phi_1 = parameters(4);

s = linspace(0,rho_1,2*nS - 1);

obj = GenImage(x,y,x0_1,y0_1,phi_1,s,sX,nS,background);
obj = obj - background;
obj = obj ./ max(obj(:)) * (max(imageFrame(:))-background);
obj = obj + background;
f = sum((imageFrame(:)-obj(:)).^2);


inform = int64(1);

end