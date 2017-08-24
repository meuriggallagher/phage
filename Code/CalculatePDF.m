function psi = CalculatePDF(Pe,T,dt,dp,n,m,alpT,alpT_th, ...
        alpP,alpP_ph,Dtt,Dtt_t,Dpp)

%% Bulk
% 1 < i < n+1, 1 < j < m
% psi_(i,j) terms
ij_entries = Pe * cot(T(:)).*alpT + Pe * alpT_th ...
    + Pe./sin(T(:)).*alpP_ph + 2*Dtt/dt^2 + 2*Dpp./(dp^2*(sin(T(:)).^2));
% remove i = 1
ij_entries(1:m) = 0;
% remove i = n
ij_entries(n*m + 1 : end) = 0;

A = spdiags(ij_entries,0,(n+1)*m,(n+1)*m);

% psi_(i,j+1) terms
ij_entries = Pe./sin(T(:))/2/dp .* alpP - 1./(dp^2 * sin(T(:)).^2).*Dpp;
% remove j = m
ij_entries((1:n+1)*m) = 0;
% shift for top diagonal
ij_entries = [0;ij_entries(1:end-1)];

A = A + spdiags(ij_entries,1,(n+1)*m,(n+1)*m);

% psi_(i,j-1) terms[zeros((n+1)*m,1) ; 1]
ij_entries = -Pe./sin(T(:))/2/dp .* alpP - 1./(dp^2 * sin(T(:)).^2).*Dpp;
% remove j = 1
ij_entries((0:n)*m+1) = 0;

A = A + spdiags(ij_entries(2:end),-1,(n+1)*m,(n+1)*m);

% psi_(i+1,j) terms
ij_entries = Pe/2/dt*alpT - cot(T(:))/2/dt.*Dtt - 1/2/dt*Dtt_t - Dtt/dt/dt;
% remove i = 1
ij_entries(1:m) = 0;
% remove i = n
ij_entries(n*m + 1 : end) = 0;
% shift for mth diagonal
ij_entries = [zeros(m,1); ij_entries];

A = A + spdiags(ij_entries,m,(n+1)*m,(n+1)*m);

% psi_(i-1,j) terms
ij_entries = - Pe/2/dt*alpT + cot(T(:))/2/dt.*Dtt ...
    + 1/2/dt*Dtt_t - Dtt/dt/dt;
% remove i = 1
ij_entries(1:m) = 0;
% remove i = n
ij_entries(n*m + 1 : end) = 0;

A = A + spdiags(ij_entries(m+1:end),-m,(n+1)*m,(n+1)*m);

%% LEFT: 1 < i < n+1, j = 1
inds = (0:n)*m + 1;
dumT = T(:);
% psi_(i,m) terms
ij_entries = -Pe./sin(dumT(inds))/2/dp .* alpP(inds) ...
    - 1./(dp^2 * sin(dumT(inds)).^2).*Dpp(inds);

A = A + sparse(inds,inds + m-1, ij_entries,(n+1)*m,(n+1)*m);

%% RIGHT: 1 < i < n+1, j = m
inds = (1:n+1)*m;
% psi_(i,1) terms
ij_entries = Pe./sin(dumT(inds))/2/dp .* alpP(inds) ...
    - 1./(dp^2 * sin(dumT(inds)).^2).*Dpp(inds);

A = A + sparse(inds,inds-(m-1),ij_entries,(n+1)*m,(n+1)*m);


%% TOP: i = 1, j = 2:m-1
inds = 2:m-1;

% psi_(1,j) terms
ij_entries = Pe * cot(dumT(inds)).*alpT(inds) + Pe * alpT_th(inds) ...
    + Pe./sin(dumT((inds))).*alpP_ph(inds) - Dtt(inds)/dt^2 ...
    + 2*Dpp(inds)./(dp^2*(sin(dumT(inds)).^2)) ...
    - Pe/dt * alpT(inds) + cot(dumT(inds))/dt.*Dtt(inds) ...
    + 1/dt*Dtt_t(inds);

A = A + sparse(inds,inds,ij_entries,(n+1)*m,(n+1)*m);

% psi_(2,j) terms
ij_entries = Pe/dt*alpT(inds) - cot(dumT(inds))/dt.*Dtt(inds) ...
    -1/dt*Dtt_t(inds) + 2/dt/dt*Dtt(inds);

A = A + sparse(inds,inds+m,ij_entries,(n+1)*m,(n+1)*m);

% psi_(3,j) terms
ij_entries = -1/dt/dt*Dtt(inds);

A = A + sparse(inds,inds+2*m,ij_entries,(n+1)*m,(n+1)*m);

%% BOTTOM i = n+1, j = 2:m-1
inds = n*m+2:(n+1)*m-1;

% psi_(n+1,j) terms
ij_entries = Pe * cot(dumT(inds)).*alpT(inds) + Pe * alpT_th(inds) ...
    + Pe./sin(dumT((inds))).*alpP_ph(inds) - Dtt(inds)/dt^2 ...
    + 2*Dpp(inds)./(dp^2*(sin(dumT(inds)).^2)) ...
    - Pe/dt * alpT(inds) + cot(dumT(inds))/dt.*Dtt(inds) ...
    + 1/dt*Dtt_t(inds);

A = A + sparse(inds,inds,ij_entries,(n+1)*m,(n+1)*m);

% psi_(n,j) terms
ij_entries = Pe/dt*alpT(inds) - cot(dumT(inds))/dt.*Dtt(inds) ...
    -1/dt*Dtt_t(inds) + 2/dt/dt*Dtt(inds);

A = A + sparse(inds,inds-m,ij_entries,(n+1)*m,(n+1)*m);

% psi_(n-1,j) terms
ij_entries = -1/dt/dt*Dtt(inds);

A = A + sparse(inds,inds-2*m,ij_entries,(n+1)*m,(n+1)*m);

%% TOP LEFT: i = 1, j = 1
% psi_(1,1)
A(1,1) = Pe*cot(T(1))*alpT(1) + Pe*alpT_th(1) + Pe/sin(T(1))*alpP_ph(1) ...
    - Dtt(1)/dt/dt - Pe/dt*alpT(1) + cot(T(1))/dt*Dtt(1) + 1/dt*Dtt_t(1)...
    + 2/dp/dp/sin(T(1))/sin(T(1))*Dpp(1);

% psi_(2,1)
A(1,1+m) = Pe/dt*alpT(1) - cot(T(1))/dt*Dtt(1) - 1/dt*Dtt_t(1) ...
    + 2/dt/dt*Dtt(1);

% psi_(3,1)
A(1,1+2*m) = -1/dt/dt*Dtt(1);

% psi_(1,2) already complete
% psi_(1,m) already complete

%% TOP RIGHT: i = 1, j = m
% psi_(1,m)
A(m,m) = Pe*cot(T(1))*alpT(1) + Pe*alpT_th(1) + Pe/sin(T(1))*alpP_ph(1) ...
    - Dtt(1)/dt/dt - Pe/dt*alpT(1) + cot(T(1))/dt*Dtt(1) + 1/dt*Dtt_t(1)...
    + 2/dp/dp/sin(T(1))/sin(T(1))*Dpp(1);

% psi_(2,m)
A(m,2*m) = Pe/dt*alpT(1) - cot(T(1))/dt*Dtt(1) - 1/dt*Dtt_t(1) ...
    + 2/dt/dt*Dtt(1);

% psi_(3,m)
A(m,3*m) = -1/dt/dt*Dtt(1);

% psi_(1,m-1) already complete
% psi_(1,1) already complete

%% BOTTOM LEFT: i = n+1, j = 1
% psi_(n+1,1)
A(m*n+1,m*n+1) = Pe*cot(T(end))*alpT(m*n+1) + Pe*alpT_th(m*n+1) ...
    + Pe/sin(T(end))*alpP_ph(m*n+1) - 1/dt/dt*Dtt(end) ...
    - Pe/dt*alpT(m*n+1) + cot(T(end))/dt*Dtt(end) + 1/dt*Dtt_t(end) ...
    + 2/dp/dp/sin(T(end))/sin(T(end))*Dpp(end);

% psi_(n,1)
A(m*n+1,m*(n-1)+1) = Pe/dt*alpT(m*n+1) + cot(T(end))/dt*Dtt(end) ...
    - 1/dt*Dtt_t(end) + 2/dt/dt*Dtt(end);

% psi_(n-1,1)
A(m*n+1,m*(n-2)+1) = -1/dt/dt*Dtt(end);

% psi_(n+1,2) already complete
% psi_(n+1,m) already complete

%% BOTTOM RIGHT: i = n+1, j = m
% psi_(n+1,m)
A(m*(n+1),m*(n+1)) = Pe*cot(T(end))*alpT(m*n+1) + Pe*alpT_th(m*n+1) ...
    + Pe/sin(T(end))*alpP_ph(m*n+1) - 1/dt/dt*Dtt(end) ...
    - Pe/dt*alpT(m*n+1) + cot(T(end))/dt*Dtt(end) + 1/dt*Dtt_t(end) ...
    + 2/dp/dp/sin(T(end))/sin(T(end))*Dpp(end);

% psi_(n,m)
A(m*(n+1),m*n) = Pe/dt*alpT(m*n+1) + cot(T(end))/dt*Dtt(end) ...
    - 1/dt*Dtt_t(end) + 2/dt/dt*Dtt(end);

% psi_(n-1,m)
A(m*(n+1),m*(n-1)) = -1/dt/dt*Dtt(end);

%% Normalisation condition
%options = optimoptions('lsqlin','display','none', ...
%    'constrainttolerance',1e-14,'optimalitytolerance',1e-14);

%psi = lsqlin(A,zeros((n+1)*m,1),[],[], ...
%    dt*dp*sin(T(:)'),1,zeros((n+1)*m,1),[],[],options);
A = [A; dt*dp*sin(T(:)')];
psi = A\[zeros((n+1)*m,1);1];
psi = reshape(psi,n+1,m);

end