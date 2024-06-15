function [pe]=PE_hassibi(p,k,n)

pe=(2^k)*exp(-(p)*dmin(k,n));

end

function [out]=dmin(k,n)
out=n*GV(1-(k/n));
end


function [out]=GV(r)
syms x
f(x)=-x*log2(x)-(1-x)*log2(1-x);

out=real(vpasolve(f(x)==r));
end

