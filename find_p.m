function [ash,sun]=find_p(n,k)

theta=.1;

p=.01:.01:.99;
pe_ash=zeros(length(p),1);
pe_sun=zeros(length(p),1);
epsilon=p.^k;
for i=1:1:length(p)
    pe_ash(i)=PE_ashikhmin(epsilon(i),k,n);
    pe_sun(i)=PE_hassibi(epsilon(i),k,n);
end

constant=theta/(n^k)*ones(length(p),1);

ash=intersections(p',pe_ash,p',constant,1);
sun=intersections(p',pe_sun,p',constant,1);

% plot(p, log(pe_sun))
% hold on
% plot(p, log(constant))
end
    
    