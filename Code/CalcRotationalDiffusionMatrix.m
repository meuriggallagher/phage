function D=CalcRotationalDiffusionMatrix(theta,phi,xbf,Xbf,X0,eps,domain,blockSize)

% Calculates the rotational diffusion matrix for an object from its body
% frame points and orientation given in spherical polars
% input:  theta, phi spherical polar angles
%         xbf        body frame collocation/force points
%         Xbf        body frame stokeslet points
%         X0         centre of rotation
%         eps        numerical regularisation length
%         domain     'i' for infinite fluid, 'h' for half space above no-slip boundary (x3>0)
%         blockSize  controls matrix assembly blocking size for numerical efficiency, 0.2 is a safe choice
% output: D          2x2 matrix of diffusion coefficients 
%                    (scaled wrt kB T / mu) in spherical polars (th,phi)

n_th=length(theta);
theta
D=zeros(2,2,n_th);

for j=1:n_th
    j,n_th
    th=theta(j);
    [e_r e_th e_phi]=BasisVectorsSphericalPolars(th,phi);
    Om_phi=e_phi;  %rotation mode in phi direction
    Om_th =e_th;   %rotation mode in theta direction
    gamma=0.0;     %no incident flow

    x=RotateToSphericalPolars(th,phi,X0,xbf);
    X=RotateToSphericalPolars(th,phi,X0,Xbf);    

    MM_phi=SolveTetheredResistanceProblem(x,X,X0,@ShearFlow,gamma,Om_phi,eps,domain,0.2);
    MM_th =SolveTetheredResistanceProblem(x,X,X0,@ShearFlow,gamma,Om_th, eps,domain,0.2);

    R(1,1)=MM_th'*e_th;
    R(1,2)=MM_phi'*e_th;
    R(2,1)=MM_th'*e_phi;
    R(2,2)=MM_phi'*e_phi;
    
    D(:,:,j)=inv(R);
end
