function [Om]=SolveTetheredMobilityProblem(x,X,X0,Uinf,Uparam,MM,eps,domain,blockSize)

  % Solves the problem of finding the angular velocity of a body constrained so that its velocity is of the form omega x (X-X0)
  % Input:  x          collocation/force grid
  %         X          stokeslet grid
  %         X0         origin of rotation
  %         Uinf       function yielding undisturbed flow velocity -- pass this function using e.g. @ShearFlow
  %         Uparam     velocity parameter e.g. (d/dt) gamma
  %         MM         imposed moment (e.g. Brownian torque) about X0
  %         eps        regularisation parameter
  %         domain     'i' infinite fluid (regularised stokeslet) or 'h' halfspace x3 > 0 (regularised blakelet)
  %         blockSize  maximum memory storage for stokeslet matrix in GB
  % Output: Om         angular velocity

  tic

  M=length(x)/3;

  % lhs assembly
  [AS NN]=AssembleStokesletMatrix(x,X,x,eps,domain,blockSize);

  [x1 x2 x3]=ExtractComponents(x);ze=0*x1;
  x1=x1-X0(1);
  x2=x2-X0(2);
  x3=x3-X0(3);

  AOm=[ze -x3 x2; x3 ze -x1; -x2 x1 ze];

  [X1 X2 X3]=ExtractComponents(X'*NN);Ze=0*X1;
  X1=X1-X0(1);
  X2=X2-X0(2);
  X3=X3-X0(3);

  AM=[Ze -X3 X2; X3 Ze -X1; -X2 X1 Ze];

  A=[AS AOm; AM zeros(3,3)];

  % rhs assembly
  [u1,u2,u3]=Uinf(x1,x2,x3,Uparam);
  rhs=[-u1;-u2;-u3;MM];

  'system assembly:'
  toc

  % solve system
  tic

  sol=A\rhs;
  f=sol(1:3*M);
  Om=sol(3*M+1:3*M+3);

  'system solver:'
  toc
