function[converged,L_b,R_b]= Recovery1(m,n,r,p_obs,iterations,error_threshold,A, obs_pattern,C0,C1,L_b,R_b,L_b_fixed_vals)

% m=100;
% n=100;
% r=5;
% p_obs=.1;
% iterations=50;
% error_threshold=10^-6; % for BP convergence

% L_i=randi(2,m,r)-1;
% R_i=randi(2,r,n)-1;
% A=mod(L_i*R_i,2);

%
% obs_pattern=rand(m,n);
% obs_pattern=(obs_pattern<p_obs);
%
%
% C0=zeros(2^r,2^r);
% for i=1:1:2^r
%     for j=1:1:2^r
%         if (mod(de2bi(i-1,r)*de2bi(j-1,r)',2)==0)
%             C0(i,j)=1;
%         end
%     end
% end
%
%
% C1=1-C0;






% Constructing the graph

% L_b=ones(2^r,m)/2^r;
% R_b=ones(2^r,n)/2^r;

%L_b(:,1)=de2bi(4,2^r)';


[~,num_of_fixed]=size(L_b_fixed_vals);
L_b(:,1:num_of_fixed)=L_b_fixed_vals;

converged=false;
for iter=1:1:iterations
    L_b_new=ones(2^r,m);
    R_b_new=ones(2^r,n);
    for i=1:1:m
        for j=1:1:n
            
            if (obs_pattern(i,j)==1)
                
                if(A(i,j)==1)
                    constraint=C1;
                else
                    constraint=C0;
                end
                
                %update L
                
                temp=constraint*R_b(:,j);
                %tempp=prod(temp','native');
                L_b_new(:,i)=temp.*L_b_new(:,i);
                L_b_new(:,i)=L_b_new(:,i)/sum(L_b_new(:,i));
                
                %update R
                
                temp=constraint*L_b(:,i);
                %tempp=prod(temp','native');
                R_b_new(:,j)=temp/sum(temp).*R_b_new(:,j);
                R_b_new(:,j)=R_b_new(:,j)/sum(R_b_new(:,j));
            end
        end
    end
    L_b1=L_b_new./sum(L_b_new);
    R_b1=R_b_new./sum(R_b_new);
    
    if  norm(L_b-L_b1)+norm(R_b-R_b1)<error_threshold
        
        converged=true;
        break
    end
    
    L_b=L_b_new./sum(L_b_new);
    L_b(isnan(L_b))=0;
    L_b(:,1:num_of_fixed)=L_b_fixed_vals;
    R_b=R_b_new./sum(R_b_new);
    R_b(isnan(R_b))=0;
end
L_b;
iter;
R_b;

end














