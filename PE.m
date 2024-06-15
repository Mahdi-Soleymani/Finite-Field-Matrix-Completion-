function [pe]=PE(epsilon,k,n)

pe=0;
for i=1:1:n
    pe=pe+nchoosek(n,i)*epsilon^i*((1-epsilon)^(n-i))*prod(n,k,i)
end
end

function [out]=prod(n,k,i)
out=1;
for j=1:1:k
    out=out*(1-2^(j-1-n+i))
end
end
