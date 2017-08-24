function alpha=CalcAdvectionUnitShear(Th,Phi,xbf,Xbf,X0,eps,domain,blockSize)

% Calculates rate of change of direction vector in spherical polars
% in unit shear flow
% under zero moment balance
% alpha = Om x e_r
% alpha . e_th  = Om_phi  =  Om . phi
% alpha . e_phi = -Om_th = -Om . th

% input:  Th, Phi      spherical polar angles, supply as column vectors of equal length
%         xbf          body frame collocation/force points
%         Xbf          body frame stokeslet points
%         X0           centre of rotation
%         eps          numerical regularisation length
%         domain       'i' for infinite fluid, 'h' for half space above no-slip boundary (x3>0)
%         blockSize    controls matrix assembly blocking size for numerical efficiency, 0.2 is a safe choice
% output: alpha(1,:,:) rate of change of direction vector, e_th component
%         alpha(2,:,:) rate of change of direction vector, e_phi component

gamma=1;  %unit shear flow
MM=[0;0;0];     %moment-free

for ii=1:size(Th,1)
    for jj=1:size(Phi,2)
    	ii
        Th(ii,jj)
        Phi(ii,jj)
        x=RotateToSphericalPolars(Th(ii,jj),Phi(ii,jj),X0,xbf);
        X=RotateToSphericalPolars(Th(ii,jj),Phi(ii,jj),X0,Xbf);

        Om=SolveTetheredMobilityProblem(x,X,X0,@ShearFlow,gamma,MM,eps,domain,blockSize)

        [e_r e_th e_phi]=BasisVectorsSphericalPolars(Th(ii,jj),Phi(ii,jj));
        alpha(1,ii,jj)=dot(e_phi,Om);
        alpha(2,ii,jj)=dot(-e_th,Om);
    end
end
