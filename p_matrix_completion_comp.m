function []=p_matrix_completion_comp()
k=2;

n=20:1:50;

ash=zeros(1,length(n));
sun=zeros(1,length(n));

for i=1:1:length(n)
    i
    [as,su]=find_p(n(i),k);
    if length(as)==0
        ash(i)=1
    else
    ash(i)=as;
    end
    
    if length(su)==0
    sun(i)=1;
    else
    sun(i)=su;
    end
end

plot(n,ash);
hold on
plot(n, sun);
end


