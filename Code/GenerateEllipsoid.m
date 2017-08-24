function X=GenerateEllipsoid(nth,nt,arad,len,centre)

th=linspace(0,2*pi,nth+1)';
th(end)=[];

tt=linspace(-len/2,len/2,nt+2);
tt(end)=[];
tt(1)=[];

X1=arad*cos(th)*sqrt(1-(2*tt/len).^2);
X2=arad*sin(th)*sqrt(1-(2*tt/len).^2);
X3=ones(nth,1)*tt;

X1=X1(:);
X2=X2(:);
X3=X3(:);

% X1=[X1;0];
% X2=[X2;0];
% X3=[X3;-len/2];
% 
% X1=[X1;0];
% X2=[X2;0];
% X3=[X3;len/2];

switch centre
    case 0
        X3=X3+len*0.5;
end

X=[X1;X2;X3];
