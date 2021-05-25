function [ b ] = construct_laguerre_bases( im_size )

%% this implements the laguerre basis as described in paper on laguerre expansion of flim

alpha =.9;

L=10; %% max order of Laguerre poly
B=[];
%%Does not use the actual data, same for all vectors
for k=1:im_size
    for l=0:L-1
        
        B(k,l+1)=(alpha).^((k-l)/2)*(1-alpha).^.5;
        sum1=0;
        for m=0:l
            
            if(m<=k)
                
                sum1=sum1 + (-1)^m*nchoosek(k,m)*nchoosek(l,m)*(alpha)^(l-m)*(1-alpha)^m;
                
            else
                sum1=sum1+ (-1)^m*nchoosek(l,m)*(alpha)^(l-m)*(1-alpha)^m;%%should this be divided by fact(m)
            end
        end
        B(k,l+1)=B(k,l+1)*sum1;
        
    end
end
       


b=reshape(B(:,:),[im_size L]);


end

