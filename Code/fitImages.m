%FITIMAGES Fits simulated bacteriophage to set of experimental images
%   FITIMAGES(FILENAME) Fits simulated bacteriophage to the set of
%   experimental images FILENAME
%
%   INPUTS:
%       FILENAME - .tif file of experimental images
%
%   OUTPUTS:
%       Output is a .mat file called FILENAME containing the processed
%       experimental images, and parameters of fitted bacteriophage
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function fitImages(filename)

% Load images
[img,nFrame,frames] = LoadImage(filename);
image = zeros(size(img{frames(1)},1),size(img{frames(1)},2),nFrms);
for ii = 1:nFrms
    dum = img{frames(ii)};
    image(:,:,ii) = double(dum);
end

%% Set up 
fprintf('Frames %i to %i...',frames(1),frames(end))

savename = [filename(1:end-4) '.mat'];

fileDetails = ['This file fits the phage and optical parameters ' ...
    'to experimental images. The aim of this fit is to find the ' ...
    'length of projected phage with theta = pi/2']; %#ok<NASGU>

parameterStrings = {'x0';'y0';'phi';'rho';'sX'; ...
    'background';'error'}; %#ok<NASGU>

paramsFit = NaN(7,nFrms);

processedImages = cell(nFrms,1);

%% Perform fit
for ii = 1:nFrame
    
    % If frame too noisy then skip
    dum = image(:,:,ii);
    if median(dum(:))/max(dum(:)) > 0.75
        paramsFit(:,ii) = inf(7,1);
        processedImages{jj} = image(:,:,ii);
        continue
    end
    
    % Preprocessing step
    imageFrame = image(:,:,ii);
    imageFrame = medfilt2(imageFrame,[5 5]);
    imageFrame = imageFrame - median(imageFrame(:));
    imageFrame(imageFrame < 3) = 0;
    
    processedImages{ii} = imageFrame;

    % Fit image
    paramsFit(:,ii) = Phage2DFit(imageFrame,remainingFrms);

end

fprintf('saving...')
save(savename,'fileDetails','parameterStrings','paramsFit', ...
    'processedImages','-v7.3')

fprintf('complete.\n')

end