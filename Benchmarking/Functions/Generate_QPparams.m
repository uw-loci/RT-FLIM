function x = Generate_QPparams(H,f,A,b,G,h,m,Label)
%%% sum(sum((V*X-B).*(V*X - B))) == trace(X'*V'*V*X-2*B'*V*X+B'*B)


%%M is the number of columns of the variable matrix. 
if(nargin<1)
%%variable should be n*m
    n=10;
    m=64;
    
    %%reading presaved mat file
    
%      load H.mat
%     H=HH;
%     G=D;
%     A=[];
%     b=[];
    alpha=0;
    V=rand(n-2,n);
    [U S]=eig(V'*V);
    
    H=U*abs(S)*U';
    H=(H+H')/2;
    f=rand(m,n);
    A=zeros(1,n);
    b=zeros(1,1);
    G=rand(4,n);
    G=ones(1,4)*G;
    h=zeros(1,1);
    Label=sign(rand(m,1)-.5);
    Label(find(Label==-1))=0;
    L=(Label*Label');
    LG=diag(sum(L,1))-L;
    Label=reshape(Label,[8 8]);
    %dlmwrite('label.txt',Label,',');
    %H=H+alpha*LG; %objective is tr(alpha'*H*alpha) + trace(alpha*L*alpha')
    %where alpha is what we are solving for
end
if(nargin>4)
    G=ones(1,size(G,1))*G;
end
if(m>1)
    n=size(H,1);
    
    %Hkron = kron(eye(m),H);
    f1=f';
    f_save=f;
    f=f1(:);
    A=zeros(1,n*m);
    b=zeros(1,1);
    
    G_save=G;
    G=(kron(eye(m),G));
    
    
   
end
% [size(H) size(f)]
% pause
h=zeros(size(G,1),1);
    
n=size(H,1);

if(numel(A)==0)
    A=zeros(1,n);
end
if(numel(b)==0)
    b=zeros(1,1);
end

A(end)=.0000000001;
HH=H;
if(m==1)
    D=G;
else
    D=G';
    
end

x=rand(n,m);

end