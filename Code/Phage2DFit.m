%PHAGE2DFIT Fits a bacteriophage to experimental images
%   FITPARAMS = PHAGE2DFIT(IMAGE) Calculates the parameters FITPARAMS which
%   best fit the IMAGE data
%
%   INPUTS:
%       IMAGE - Experimental image data
%
%   OUTPUTS:
%       FITPARAMS - [x0;y0;phi;rho;sX;background;error]
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function fitParams = Phage2DFit(image)

global imageFrame x y x0 y0 phi rho sX background nS

%% Image parameters
[nY,nX] = size(image);
[~,~,~,x,y] = ImageParameters(nX,nY);
nS = 35; % Phage discretisation
fitParams = nan(7,1);
imageFrame = image;
background = 0;

%% Nag options
nagParams = NagGlobSetup;

%% First initialise x0, y0,rho bounds based on high intensity part of image
[~,dum] = max(image(:));dx = x(:); dx = dx(dum);dy = y(:); dy = dy(dum);

% Upper and lower bounds for phage parameters
%           x0     y0      rho  phi
blPhage = [dx-1e-6,dy-1e-6,1e-9,-0.5];
buPhage = [dx+1e-6,dy+1e-6,1e-6,0.5+pi];


%% Fit x0, y0, phi, rho with a guess of sX
sX = 4*ImagingParameters;

[~, ~, ~, ~, ~, params,~,~,~,iFail] = nag_glopt_bnd_mcs_solve( ...
    'MinFuncPhage',nagParams{1,3}, blPhage, buPhage, ...
    nagParams{4,1}, nagParams{4,2},nagParams{4,3}, nagParams{1,1}, ...
    nagParams{1,2});

% Alter boundary conditions if fit fails
runcount = 1;
while iFail == 7 || ...
        (min(abs((blPhage - params)/blPhage)) < 1e-1) || ...
        (min(abs((buPhage - params)/buPhage)) < 1e-1)
    blPhage(4) = blPhage(4) + 0.5;
    buPhage(4) = buPhage(4) + 0.5;
    [~, ~, ~, ~, ~, params,~,~,~,iFail] = nag_glopt_bnd_mcs_solve( ...
    'MinFuncPhage',nagParams{1,3}, blPhage, buPhage, ...
    nagParams{4,1}, nagParams{4,2},nagParams{4,3}, nagParams{1,1}, ...
    nagParams{1,2});

    runcount = runcount + 1;
    if runcount == 8
        return
    end
end
    
if iFail ~= 0
    return
end

%% Fit Gaussian point spread function width sX
x0 = params(1);
y0 = params(2);
rho = params(3);
phi = params(4);

%            sX
blOptics = sX/10;
buOptics = sX*10;

% Solve
[~, ~, ~, ~, ~, params,~,~,~,iFail] = nag_glopt_bnd_mcs_solve( ...
    'MinFuncOpts',nagParams{1,3}, blOptics, buOptics, nagParams{11,1},...
    nagParams{11,2},nagParams{11,3}, nagParams{1,1}, nagParams{1,2});

if iFail ~= 0
    return
end

%% Refit x0 y0 phi theta with new sX
sX = params(1);

[~, ~, ~, ~, ~, params,~,~,~,iFail] = nag_glopt_bnd_mcs_solve( ...
    'MinFuncPhage',nagParams{1,3}, blPhage, buPhage, nagParams{4,1},...
    nagParams{4,2},nagParams{4,3}, nagParams{1,1}, nagParams{1,2});

% Alter boundary conditions if fit fails
runcount = 1;
while iFail == 7 || ...
        (min(abs((blPhage - params)/blPhage)) < 1e-1) || ...
        (min(abs((buPhage - params)/buPhage)) < 1e-1)
    blPhage(4) = blPhage(4) + 0.5;
    buPhage(4) = buPhage(4) + 0.5;
    [~, ~, ~, ~, ~, params,~,~,~,iFail] = nag_glopt_bnd_mcs_solve( ...
    'MinFuncPhage',nagParams{1,3}, blPhage, buPhage, ...
    nagParams{4,1}, nagParams{4,2},nagParams{4,3}, nagParams{1,1}, ...
    nagParams{1,2});

    runcount = runcount + 1;
    if runcount == 8
        return
    end
end

if iFail ~= 0
    return
end

%% Refit sX background
x0 = params(1);
y0 = params(2);
rho = params(3);
phi = params(4);

%             sX     background
blIntsopts = [sX/4 ; 0];
buIntsopts = [sX*4;  max(imageFrame(:))/4];

% Solve
[~, ~, ~, ~, ~, params,~,~,~,iFail] = nag_glopt_bnd_mcs_solve( ...
    'MinFuncIntsOpts',nagParams{1,3}, blIntsopts, buIntsopts, ...
    nagParams{2,1},nagParams{2,2},nagParams{2,3}, nagParams{1,1}, ...
    nagParams{1,2});

if iFail ~= 0
    return
end

%% Refit x0 y0 phi rho with new sX, background
sX = params(1);
background = params(2);

[~, ~, ~, ~, ~, params,error,~,~,iFail] = nag_glopt_bnd_mcs_solve( ...
    'MinFuncPhage',nagParams{1,3}, blPhage, buPhage, ...
    nagParams{4,1}, nagParams{4,2},nagParams{4,3}, nagParams{1,1}, ...
    nagParams{1,2});

% Alter boundary conditions if fit fail
runcount = 1;
while iFail == 7 || ...
        (min(abs((blPhage - params)/blPhage)) < 1e-1) || ...
        (min(abs((buPhage - params)/buPhage)) < 1e-1)
    blPhage(4) = blPhage(4) + 0.5;
    buPhage(4) = buPhage(4) + 0.5;
    [~, ~, ~, ~, ~, params,~,~,~,iFail] = nag_glopt_bnd_mcs_solve( ...
    'MinFuncPhage',nagParams{1,3}, blPhage, buPhage, ...
    nagParams{4,1}, nagParams{4,2},nagParams{4,3}, nagParams{1,1}, ...
    nagParams{1,2});

    runcount = runcount + 1;
    if runcount == 8
        return
    end
end

if iFail ~= 0
    return
end

%% Save
fitParams(1) = params(1); % x0
fitParams(2) = params(2); % y0
fitParams(3) = params(4); % phi
fitParams(4) = params(3); % projected length
fitParams(5) = sX;
fitParams(6) = background;
fitParams(7)= error;

end