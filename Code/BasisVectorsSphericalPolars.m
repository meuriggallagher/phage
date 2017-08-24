function [e_r e_th e_phi]=BasisVectorsSphericalPolars(th,phi)

e_r=[cos(phi)*sin(th);sin(phi)*sin(th);cos(th)];
e_phi=[-sin(phi);cos(phi);0];
e_th =[cos(phi)*cos(th);sin(phi)*cos(th);-sin(th)];