%% MTG Phage
%% Generate rod
len = 1.0; % phage length
arad = 7/900; % phage radius
aenv = 80.0; %envelope 
centre = 0; 
shape = 'rod';

nth = 4; % number of body frame collocation/force points
Nth = 9; % number of body frame stokeslet points
axisFactor = 12; %??

nt=nth*axisFactor;
Nt=Nth*axisFactor;

% force grid
xbf=GenerateRod(nth,nt,arad,aenv,len,centre);
% quadrature grid
Xbf=GenerateRod(Nth,Nt,arad,aenv,len,centre);

switch centre
    case 0
        X0=[0;0;arad]; %centre of rotation shifted avoid boundary collision
    case 1
        X0=[0;0;0];    %centre of rotation is the origin
end

eps=0.01; % small parameter for reg. stokeslets
domain='h'; % hemispherical domain

% controls matrix assembly blocking size for numerical efficiency, 0.2 safe
blockSize=0.2; 

switch domain
    case 'h'
        thmax = pi/2 * 0.98;
    case 'i'
        thmax = pi * 0.99;
end

n = length(xbf)/3;
N = length(Xbf)/3;

%% Calculate alpha (rotational advection vector)
Nth=10;
Nphi=10;

Nth1=100;
Nphi1=100;

[ppAlpha,~,phigrid,alphagrid,thgrid1,phigrid1,alphagrid1] = ...
    CalcSplineAdvectionUnitShear_MTG(Nth,Nphi,Nth1,Nphi1,thmax,xbf,Xbf, ...
    X0,eps,domain,blockSize);

[Phigrid1,Thgrid1]=meshgrid(phigrid1(:),thgrid1(:));

%% Calculate D (rotational diffusion matrix)
[ppD,thgrid,Dgrid]=CalcSplineRotationalDiffusionMatrix(Nth,thmax,xbf,...
    Xbf,X0,eps,domain,blockSize);

%% Calculate PDFs
% Set up
n = 99; m = 100;
th = linspace(0,pi/2,n+3); th = th(2:end-1);dt = th(2)-th(1);
ph = linspace(-pi,pi,m+1); ph = ph(1:end-1); dp = ph(2)-ph(1);
[T,P] = meshgrid(th,ph);


% calculate advection and diffusion (and derivatives)
alpT = fnval(ppAlpha{1},[T(:)';P(:)'])';
alpT_th = fnval(fnder(ppAlpha{1},[1;0]),[T(:)';P(:)'])';
alpT_ph = fnval(fnder(ppAlpha{1},[0;1]),[T(:)';P(:)'])';

alpP = fnval(ppAlpha{2},[T(:)';P(:)'])';
alpP_th = fnval(fnder(ppAlpha{2},[1;0]),[T(:)';P(:)'])';
alpP_ph = fnval(fnder(ppAlpha{2},[0;1]),[T(:)';P(:)'])';

Dtt = fnval(ppD(1,1),T(:)')';
Dtt_t = fnval(fnder(ppD(1,1)),T(:)')';
Dpp = fnval(ppD(2,2),T(:)')';

% Peclet number
Pe = linspace(1e-2,500,1e3);

% PDF
psi = cell(length(Pe),1);

% Calculate PDFs
for ii = 1:length(Pe)
    if mod(ii,100) == 0 
        fprintf('Solving Pe = %e\n',Pe(ii)) 
    end
    psi{ii} = CalculatePDF(Pe(ii),T,dt,dp,n,m,alpT,alpT_th, ...
        alpP,alpP_ph,Dtt,Dtt_t,Dpp);
end

%% Calculate marginal PDF
% Set up points for integration
dt = (T(1,2)-T(1,1))/2;
dp = (P(2,1)-P(1,1))/2;
t = min(T(:)) : dt : max(T(:));
p = min(P(:)) : dp : max(P(:)); p = p(1:end-1);

% Memory allocation
dumPsi = zeros(length(p),length(inds));
[t,p] = meshgrid(t,p);

% For each Pe integrate in theta to obtain marginal PDF for phi
for ii = 1:length(Pe)
    dum = psi{ii}; % Calculated PDF
    
    % Spline and sample pdf for integration
    dum = csape({T(1,:);P(:,1)'},dum','periodic');
    dum = fnval(dum,[t(:)';p(:)']);
    dum(dum < 0) = 0;
    dum = reshape(dum,size(t));
    dum = dum .* sin(t); % For integrating wrt theta
    
    % Integrate using Simpson's rule 
    d1 = dum(:,1);
    for jj = 2:size(p,1)-1
        if mod(jj,2) == 0
            d1 = d1 + 4 * dum(:,jj);
        else
            d1 = d1 + 2 * dum(:,jj);
        end
    end
    d1 = dt/3 * (d1 + dum(:,end));       
   
    % Normalise area to 1 
    sF = d1(1);
    for jj = 2:length(d1)-1
        if mod(jj,2) == 0
            sF = sF + 4*d1(jj);
        else
            sF = sF + 2*d1(jj);
        end
    end
    sF = dp/3 * (sF + d1(end));
    
    dumPsi(:,ii) = d1 / sF;
    
end

%% Spline marginal PDF in phi and Pe
marginalPDF = csapi({p(:,1)';Pe},dumPsi);

%% Calculate marginal CDF
nS = 150; % Number of points in phi
phi = linspace(-pi/2,pi/2,2*nS-1);
dP = phi(3)-phi(1);

marginalCDF = zeros(nS,length(Pe));

for ii = 1:length(Pe)
    pec = Pe(ii);
    
    % Sample PDF
    PDF = fnval(marginalPDF,[phi;pec*ones(1,length(phi))]);
    PDF(PDF < 0) = 1e-20; % Force negative values to be small
    
    % Integrate to find CDF
    for jj = 1:nS-1
        marginalCDF(jj+1,ii) = dP/6 * (PDF(2*jj-1) + ...
            4*PDF(2*jj) + PDF(2*jj+1));
    end
    marginalCDF(:,ii) = cumsum(marginalCDF(:,ii));
    
    % Scale to 1 at end
    marginalCDF(:,ii) = marginalCDF(:,ii) * 1/marginalCDF(end,ii);
    
end


marginalCDF = csapi({phi(1:2:end);Pe},marginalCDF);

save('marginalPDF.mat','marginalPDF','marginalCDF','Pe','phi','-v7.3')