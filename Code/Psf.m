%PSF Generate a Gaussian point spread function
%   PSF = PSF(X,Y,XP,YP,SX) generates a Gaussian point spread function 
%
%   INPUTS:
%       (X,Y)   - Image coordinates
%       (XP,YP) - Phage point
%       SX      - PSF width
%
%   OUTPUTS:    
%       PSF     - Gaussian point spread function
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function psf = Psf(x,y,xP,yP,sX)

psf = exp(- (x - xP).^2/(2 * sX^2) - (y - yP).^2/(2 * sX^2));

end