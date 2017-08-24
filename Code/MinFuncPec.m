%MINFUNCPEC Objective function for fitting Peclet number
%   [F,USER,INFORM]=MINFUNCPEC(N,PARAMETERS,NSTATE,USER Objective 
%   function for fitting Peclet number, set up for the NAG global 
%   optimisation regime e05jbc
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function [f,user,inform] = MinFuncPec(n,parameters,nstate,user) ...
    %#ok<INUSL>

global sampCDF marginalCDF phi

Pe = repmat(parameters,1,length(phi));

dum = fnval(marginalCDF, [phi;Pe]);

dum = abs(dum - sampCDF);

f = max(dum);

inform = int64(1);

end