function [pp,th,phi,alpha,th1,phi1,alpha1]=CalcSplineAdvectionUnitShear_MTG(Nth,Nphi,Nth1,Nphi1,thmax,xbf,Xbf,X0,eps,domain,blockSize)

% calculates cubic interpolating spline for the
% rotational diffusion matrix of a specified body

% input:  Nth        number of theta values to sample, from 0 to thmax
%         Nphi       number of phi values to sample, from 0 to 2 pi
%         thmax      maximum theta value, should be slightly less than pi/2
%         xbf        body frame collocation/force points
%         Xbf        body frame stokeslet points
%         X0         centre of rotation
%         eps        numerical regularisation length
%         domain     'i' for infinite fluid, 'h' for half space above no-slip boundary (x3>0)
%         blockSize  controls matrix assembly blocking size for numerical efficiency, 0.2 is a safe choice
% output: pp         matrix of pp-forms for the cubic spline approximation to the diffusion matrix
%         th         theta values for sampling advection
%         phi        phi values for sampling advection
%         alpha      2xNth array of sampled advection

fn_alpha=@(X,Phi) CalcAdvectionUnitShear(acos(X),Phi,xbf,Xbf,X0,eps,domain,blockSize);

switch domain
    case 'h'
        th=[linspace(0,thmax,Nth) pi/2];
        phi=linspace(0,2*pi,Nphi);

        th1=[linspace(0,thmax,Nth1) pi/2];
        phi1=linspace(0,2*pi,Nphi1);
        
    case 'i'
        th=[linspace(0,thmax,Nth) pi];
        phi=linspace(0,2*pi,Nphi);

        th1=[linspace(0,thmax,Nth1) pi];
        phi1=linspace(0,2*pi,Nphi1);
end

[Phi, Th]=meshgrid(phi,th);
[Phi1, Th1]=meshgrid(phi1,th1);

alpha=fn_alpha(cos(Th),Phi);
switch domain 
    case 'i'
    case 'h'
        alpha(:,Nth+1,:)=zeros(2,1,Nphi); % zero for max value of theta
end

pp1 = csape({Th(:,1),Phi(1,:)},squeeze(alpha(1,:,:)),'periodic');
pp2 = csape({Th(:,1),Phi(1,:)},squeeze(alpha(2,:,:)),'periodic');

Alpha1 = fnval(pp1,[Th1(:)';Phi1(:)']); Alpha1 = reshape(Alpha1,size(Th1));
Alpha2 = fnval(pp2,[Th1(:)';Phi1(:)']); Alpha2 = reshape(Alpha2,size(Th1));

pp{1}=pp1;
pp{2}=pp2;

alpha1{1}=Alpha1;
alpha1{2}=Alpha2;

