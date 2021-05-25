function [ c ] = regress_classify( b, I )

I = sum(I,3)';

n=size(b,2);
m=size(I,1);
Ii=I'; 
D=diff(b,3,1);

H=b'*b; %*2 if using quadprog


f=-2*Ii'*b;

c=Generate_QPparams(H, f, [],[], D, zeros(size(D,1), size(f,1)), m);
c=reshape(c,[n m]);