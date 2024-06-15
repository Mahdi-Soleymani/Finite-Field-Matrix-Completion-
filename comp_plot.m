function []=comp_plot()

n=50;
k=5;
p=.1:.005:.99;

hassibi=zeros(1,length(p));
ashikhmin=zeros(1,length(p));


for i=1:1:length(p)
    hassibi(i)=min(log(PE_hassibi(p(i),k,n)),0);
    ashikhmin(i)=min(log(PE_ashikhmin(p(i),k,n)),0);
%     hassibi(i)=min((PE_hassibi(p(i),k,n)),1);
%     ashikhmin(i)=min((PE_ashikhmin(p(i),k,n)),1);
end

plot(p, hassibi);
hold on
plot(p, ashikhmin);
