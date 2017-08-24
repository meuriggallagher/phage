function [u1,u2,u3]=ShearFlow(x1,x2,x3,gamma)

  M=size(x1,1);

  u1=gamma*x3;
  u2=zeros(M,1);
  u3=zeros(M,1);

