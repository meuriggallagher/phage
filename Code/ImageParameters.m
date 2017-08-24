%IMAGEPARAMETERS Outputs the image parameters
%   [PIXSIZE,LX,LY,X,Y] = IMAGEPARAMETERS(NX,NY) Outputs the image 
%   parameters for the experimental setup
%
%   INPUTS:
%       NX         - Width of image (pixels)
%       NY         - Height of image (pixels)
%
%   OUTPUTS:    
%       PIXSIZE - Size of square pixel (um)
%       LX         - Width of image (microns)
%       LY         - Height of image (microns)
%       (X,Y)      - Image coordinates
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function [pixSize,lX,lY,x,y] = ImageParameters(nX,nY)

% Pixel size
pixSize = 1e-6/14.42;

% Image size in microns
lX = nX * pixSize;
lY = nY * pixSize;

% Image coordinates
x = linspace(-lX/2,lX/2,nX);
y = linspace(-lY/2,lY/2,nY);
[x,y] = meshgrid(x,y);

end
