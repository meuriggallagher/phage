%GENPHAGE Generate a bacteriophage
%   [PX,PY] = GENPHAGE(X0,Y0,PHI,RHO) generates the location of the
%   projection of a bacteriophage in 2D space
%
%   INPUTS:
%       (X0,Y0) - Phage tethered location
%       PHI     - Phage angle
%       RHO     - Phage length
%
%   OUTPUTS:    
%       [PX,PY] - Phage location in 2D space
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function [pX,pY] = GenPhage(x0,y0,phi,rho)

pX = x0 + rho * cos(phi);

pY = y0 + rho * sin(phi);

end