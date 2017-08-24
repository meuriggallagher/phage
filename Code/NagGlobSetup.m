%NAGGLOBSETUP Setup for the NAG global optimisation regime e05jbc
%
% Copyright (c) M.T.Gallagher 2017, all rights reserved
% E-mail: m.t.gallagher@bham.ac.uk
% URL:    http://www.meuriggallagher.com/
% GIT:    https://github.com/meuriggallagher/phage
function nagParams = NagGlobSetup

nag_issue_warnings(true);

nagParams{1,1} = 'e05jbk'; % no monitoring
[nagParams{1,2}, ~] = e05ja;
nagParams{1,2} = e05jd('Function Evaluations Limit = 100000', nagParams{1,2});
nagParams{1,3} = int64(0); % All bounds will be given

for ii = 2:10
    nagParams{ii,1} = zeros(ii,3); % Only need to _declare_ the init.-list
    nagParams{ii,2} = zeros(ii,1,'int64');% these will be _set_ internally.
    nagParams{ii,3} = zeros(ii,1,'int64');
end

nagParams{11,1} = zeros(1,3); % Only need to _declare_ the init.-list
nagParams{11,2} = zeros(1,1,'int64');% these will be _set_ internally.
nagParams{11,3} = zeros(1,1,'int64');
end