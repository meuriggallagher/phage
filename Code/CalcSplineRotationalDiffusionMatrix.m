function [pp,th,D]=CalcSplineRotationalDiffusionMatrix(Nth,thmax,xbf,Xbf,X0,eps,domain,blockSize)

% calculates cubic interpolating spline for the
% rotational diffusion matrix of a specified body

% input:  Nth        number of theta values to sample, from 0 to thmax
%         thmax      maximum theta value, should be slightly less than pi/2
%         xbf        body frame collocation/force points
%         Xbf        body frame stokeslet points
%         X0         centre of rotation
%         eps        numerical regularisation length
%         domain     'i' for infinite fluid, 'h' for half space above no-slip boundary (x3>0)
%         blockSize  controls matrix assembly blocking size for numerical efficiency, 0.2 is a safe choice
% output: pp         matrix of pp-forms for the cubic spline approximation to the diffusion matrix
%         th         theta values for sampling diffusion matrices
%         D          2x2xNth array of sampled diffusion matrices

phi=0;
fn_D=@(x) CalcRotationalDiffusionMatrix(acos(x),phi,xbf,Xbf,X0,eps,domain,blockSize);

switch domain
    case 'i'
        th = [linspace(0,thmax,Nth), pi];
    case 'h'
        th=[linspace(0,thmax,Nth), pi/2];
end

%size(fn_D(cos(th)))

D=fn_D(cos(th));
switch domain
    case 'h'
        D(:,:,Nth+1)=zeros(2) + eps;
end

[pp11]=interp1(th,squeeze(D(1,1,:)),'pchip','pp');
[pp12]=interp1(th,squeeze(D(1,2,:)),'pchip','pp');
[pp21]=interp1(th,squeeze(D(2,1,:)),'pchip','pp');
[pp22]=interp1(th,squeeze(D(2,2,:)),'pchip','pp');

pp=[pp11 pp12; pp21 pp22];

