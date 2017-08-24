%GENIMAGE Simulate a fluorescent microscopy image of bacteriophage
%   IMG = GENIMAGE(X,Y,X0,Y0,PHI,S,SX,NS,BACKGROUND) generates a simulated
%   fluroescent microscopy image of a bacteriophage through convolution
%   with a Gaussian point spread function.
%
%   INPUTS:
%       X,Y        - Image coordinates
%       (X0,Y0)    - Phage tethered location
%       PHI        - Phage angle
%       S          - Phage length discretisation for Simpson's rule
%       integration
%       SX         - Gaussian PSF width
%       NS         - Number of integrals to compute
%       BACKGROUND - Either single valued background or specified ver
%       size(IMG)
%
%   OUTPUTS:
%       IMG        - Simulated image
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function img = GenImage(x,y,x0,y0,phi,s,sX,nS,background)

% Generate projected bacteriophage
[xP,yP] = GenPhage(x0,y0,phi,s);

% Step-size in arclength
ds = s(2)-s(1);

% PSF at first three points of bacteriophage
p1 = Psf(x,y,xP(1),yP(1),sX);
p2 = Psf(x,y,xP(2),yP(2),sX);
p3 = Psf(x,y,xP(3),yP(3),sX);

% Integrate using Simpson's rule
img = ds/6 * (p1 + 4*p2 + p3);

% Iterate along bacteriophage length
for ii = 2:nS-1 
    p1 = p3;
    p2 = Psf(x,y,xP(2*ii),yP(2*ii),sX);
    p3 = Psf(x,y,xP(2*ii+1),yP(2*ii+1),sX);
    
    img = img + ds/6 * (p1 + 4*p2 + p3);
end

% Add background
img = img + background;

end