function [MM]=SolveTetheredResistanceProblem(x,X,X0,Uinf,Uparam,Om,eps,domain,blockSize)

  % Solves the problem of finding the angular velocity of a body constrained so that its velocity is of the form omega x (X-X0)
  % Input:  x          collocation/force grid
  %         X          stokeslet grid
  %         X0         origin of rotation
  %         Uinf       function yielding undisturbed flow velocity -- pass this function using e.g. @ShearFlow
  %         Uparam     velocity parameter e.g. (d/dt) gamma
  %         Om         imposed angular velocity
  %         eps        regularisation parameter
  %         domain     'i' infinite fluid (regularised stokeslet) or 'h' halfspace x3 > 0 (regularised blakelet)
  %         blockSize  maximum memory storage for stokeslet matrix in GB
  % Output: MM         moment about X0 due to rotation and incident flow

  tic

  % lhs assembly
  [AS, NN]=AssembleStokesletMatrix(x,X,x,eps,domain,blockSize);

  % rhs assembly
  
  [x1, x2, x3]=ExtractComponents(x);%ze=0*x1;
  [u1,u2,u3]=Uinf(x1,x2,x3,Uparam);
  x1=x1-X0(1);
  x2=x2-X0(2);
  x3=x3-X0(3);
  rhs=[-u1;-u2;-u3]+[Om(2)*x3-Om(3)*x2;Om(3)*x1-Om(1)*x3;Om(1)*x2-Om(2)*x1]; % -u + omega cross x
  
  t = toc;
  fprintf('system assembly: time taken is %e seconds\n',t)

  % solve system
  tic
  f=AS\rhs;  
  cond(AS)
  t = toc;
  fprintf('system solver: time taken is %e seoncds\n',t)
  
  % calculate moment
  [X1, X2, X3]=ExtractComponents(X'*NN);Ze=0*X1;
  X1=X1-X0(1);
  X2=X2-X0(2);
  X3=X3-X0(3);
  AM=[Ze -X3 X2; X3 Ze -X1; -X2 X1 Ze];
  MM=AM*f;
