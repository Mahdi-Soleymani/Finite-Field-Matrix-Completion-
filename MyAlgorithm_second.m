%function [L_b,R_b]=MyAlgorithm()
function [error]=MyAlgorithm_second()



m=100;
n=100;
r=5;
p_obs=.3;
iterations=50;
error_threshold=10^-6; % for BP convergence

%Generating a random matrix
L_i=randi(2,m,r)-1;
R_i=randi(2,r,n)-1;
A=mod(L_i*R_i,2);

%Random observation pattern
obs_pattern=rand(m,n);
obs_pattern=(obs_pattern<p_obs);

L_connections=cell(1,m);
for i=1:m
    temp=[];
    for j=1:n
        if obs_pattern(i,j)==1
            temp=[temp,j];
        end
    end
    L_connections{i}=temp;
end

    

R_connections=cell(1,n);
for j=1:n
    temp=[];
    for i=1:m
        if obs_pattern(i,j)==1
            temp=[temp,i];
        end
    end
    R_connections{j}=temp;
end


%Constructing constraint matrices
C0=zeros(2^r,2^r);
for i=1:1:2^r 
    for j=1:1:2^r
        if (mod(de2bi(i-1,r)*de2bi(j-1,r)',2)==0)
            C0(i,j)=1;
        end
    end
end


C1=1-C0;


%%Initialization 



% %%% Initializing near the ground truth
% L_b_initial=zeros(2^r,m)/2^r;
% R_b_initial=zeros(2^r,n)/2^r;
% 
% for i=1:1:m
%     L_b_initial(:,i)=vector_to_cat(L_i(i,:),r);
% end
% 
% for i=1:1:n
%     R_b_initial(:,i)=vector_to_cat(R_i(:,i)',r);
% end




L_b_initial=ones(2^r,m)/2^r;

R_b_initial=ones(2^r,n)/2^r;

% L_b=ones(2^r,m)/2^r;
% 
% R_b=ones(2^r,n)/2^r;

% L_b_initial=rand(2^r,m);
% L_b_initial=L_b_initial./sum(L_b_initial);
% 
% R_b_initial=rand(2^r,n);
% R_b_initial=R_b_initial./sum(R_b_initial);

L_b_fixed_vals=[];
convergence_status=[];

for i=1:1:m
    i
%     L_b_initial=rand(2^r,m);
%     L_b_initial=L_b_initial./sum(L_b_initial)
% 
%     R_b_initial=rand(2^r,n);
%     R_b_initial=R_b_initial./sum(R_b_initial);
    
    [converged,L_b,R_b]= Recovery(m,n,r,p_obs,iterations,error_threshold,A, obs_pattern,L_connections,R_connections,C0,C1,L_b_initial,R_b_initial,L_b_fixed_vals);
    convergence_status=[convergence_status,converged]
    L_b;

    fixed_value=discretesample(L_b(:,i), 1);
    L_b_fixed_vals=[L_b_fixed_vals,double(1:2^r == fixed_value)'];
    L_b_fixed_vals;
    
    
    nnz(L_b)
    nnz(R_b)
    if nnz(L_b)<=m && nnz(R_b)<=n
        sum(L_b);
        sum(R_b);
        break
    end
    

    
end

L=zeros(m,r);


for i=1:1:m
    loc=find(L_b(:,i));
    L(i,:)=de2bi(loc(1)-1,r);
end

R=zeros(r,n);


for i=1:1:n
    loc=find(R_b(:,i));
    R(:,i)=de2bi(loc(1)-1,r)';
end
    
    
convergence_status
error=norm(A-mod(L*R,2));

    
end







function[out]=vector_to_cat(a,r)

id=bi2de(a);
out=zeros(2^r,1);
out(id+1)=1;

end







function[converged,L_b,R_b]= Recovery(m,n,r,p_obs,iterations,error_threshold,A, obs_pattern,L_connections,R_connections,C0,C1,L_b,R_b,L_b_fixed_vals)


L_b=ones(2^r,m)/2^r;
R_b=ones(2^r,n)/2^r;


[~,num_of_fixed]=size(L_b_fixed_vals);
L_b(:,1:num_of_fixed)=L_b_fixed_vals;

converged=false;
for iter=1:1:iterations
    L_b_new=ones(2^r,m);
    R_b_new=ones(2^r,n);
    
C0Rb=C0*R_b;
%C1Rb=C1*R_b;
C1Rb=1-C0Rb;
C0Lb=C0*L_b;
%C1Lb=C1*L_b;
C1Lb=1-C0Lb;
    
    
    % update L_b
    for i=1:m
        a=L_connections{i};
        
        for l=1:length(a)
            if A(i,a(l))==0
                L_b_new(:,i)=C0Rb(:,a(l)).*L_b_new(:,i);
            else
                L_b_new(:,i)=C1Rb(:,a(l)).*L_b_new(:,i);
            end
        L_b_new(:,i)=L_b_new(:,i)/sum(L_b_new(:,i));
  
        end
                
    end
    
    
    
    
    % update R_b
    for j=1:n
        a=R_connections{j};
        
        for l=1:length(a)
            if A(a(l),j)==0
                R_b_new(:,j)=C0Lb(:,a(l)).*R_b_new(:,j);
                %temp=[temp,C0Lb(:,a(l))];
            else
                R_b_new(:,j)=C1Lb(:,a(l)).*R_b_new(:,j);
                %temp=[temp,C1Lb(:,a(l))];
            end
        end
        R_b_new(:,j)=R_b_new(:,j)/sum(R_b_new(:,j));
        
    end
            
        
    
    L_b1=L_b_new./sum(L_b_new);
    R_b1=R_b_new./sum(R_b_new);
    
    if  norm(L_b-L_b1)+norm(R_b-R_b1)<error_threshold
        
        converged=true;
        break
    end
    
    L_b=L_b1;
    L_b(isnan(L_b))=0;
    L_b(:,1:num_of_fixed)=L_b_fixed_vals;
    R_b=R_b1;
    R_b(isnan(R_b))=0;
end
L_b;
iter;
R_b;

end

















    
    
    
    
    