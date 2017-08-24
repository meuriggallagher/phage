%LOADIMAGE Loads data from .tif file
%   [IMAGE,NFRAME,FRAMES] = LOADIMAGE(FILENAME) loads experimental image
%   data from .tif file for processing 
%
%   [IMAGE,NFRAME,FRAMES] = LOADIMAGE(FILENAME,FRAME) loads frames
%   [FRAME(1),FRAME(2)] from experimental image data from .tif file for 
%   processing
%
%   INPUTS:
%       FILENAME - Filename of tif
%       FRAME    - [F1,F2] frame interval to load (OPTIONAL)
%
%   OUTPUTS:
%       IMAGE    - Images loaded from .tif
%       NFRAME   - Number of frames
%       FRAMES   - Indices of frames
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function [image,nFrame,frames] = LoadImage(filename,frame)

% Image information
tiffInfo = imfinfo(filename);
nFrame = numel(tiffInfo);

% If no frames specified, take all frames
if nargin == 1
    frame = [1,nFrame];
else
    % Check first and last frame are in .tif file
    if frame(1) < 1
        frame(1) = 1;
        warning('First frame < 1, frame(1) = 1')
    end

    if frame(2) > nFrame
        frame(2) = nFrame;
        warning('End frame > nFrames, frame(2) = nFrame')
    end    
end

image = cell(nFrame,1);

for iFrame = frame(1):frame(end)
    image{iFrame} = double(imread(filename,'Index',iFrame,'Info',tiffInfo));
end

frames = frame(1):frame(end);

end