function [pe]=PE_ashikhmin(epsilon,k,n)

pe=0;
for i=1:1:n
    pe=pe+(nchoosek(n,i)*epsilon^i*((1-epsilon)^(n-i))*ash(n,k,i));
    
end
end


function [out]=ash(n,k,i)
out=0;
for r=min(0,k-n+i):1:min(k,i)
    out=out+((g(i,r)*g(n-i,k-r))/(g(n,k)))*2^(r*(n-i-k+r))*(1-2^(r-k));
end
end

function [out]=g(a,r)
out=1;
for j=0:1:r-1
    out=out*((2^a-2^j)/(2^r-2^j));
end
end