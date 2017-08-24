function X=GenerateRod(nth,nt,arad,aenv,len,centre)

%Envelope
E=@(t) 1./(1+exp(-aenv*(t)+2*len))./(1+exp(-aenv*(len-t)+2*len));

th=linspace(0,2*pi,nth+1)';
th(end)=[];

tt=linspace(0,len,nt+2);
tt(end)=[];
tt(1)=[];

X1=arad*cos(th)*E(tt);
X2=arad*sin(th)*E(tt);
X3=ones(nth,1)*tt;

X1=X1(:);
X2=X2(:);

switch centre
    case 0
        X3=X3(:);
    case 1
        X3=X3(:)-len*0.5;
end

X=[X1;X2;X3];
