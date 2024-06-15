%function [L_b,R_b]=MyAlgorithm()

function []=feeded_by_python_boolean()


recovered=readNPY('recovered.npy');
original=readNPY('original.npy');
L_mats=readNPY('L_mats.npy');
R_mats=readNPY('R_mats.npy');
masks=readNPY('masks.npy');
baseline_errors=readNPY('errors.npy');

%m,n,num_of_samples=size(masks);
num_of_samples=1;

my_errors=ones(num_of_samples,1);


m=100;
n=100;
r=5;



C0=zeros(2^r,2^r);
for i=1:1:2^r
    for j=1:1:2^r
        res=de2bi(i-1,r)*de2bi(j-1,r)'>0;
        if (res==0)
            C0(i,j)=1;
        end
    end
end


C1=1-C0;
repeat=10;
for sample=1:num_of_samples
    
for rep=1:repeat   
    rep
    
    
    
    %%% Initializing by the output of Ravanbakhsh's paper
    
    L_i=L_mats(:,:,sample);
    R_i=R_mats(:,:,sample);
    
    L_b=zeros(2^r,m)/2^r;
    R_b=zeros(2^r,n)/2^r;
    
    
    perturbation_prob=1;
    
    
    for i=1:1:m
        L_b(:,i)=vector_to_cat(L_i(i,:),r);
        %%% perturbing from the solution
        if rand<perturbation_prob
            %%% slight perturnation
            
                     L_b(:,i)=L_b(:,i)+1/2^r;
                     L_b(:,i)=L_b(:,i)/sum(L_b(:,i));

            
            %%%% replacing with random vectors
            
            %         L_b_initial(:,i)=rand(2^r,1);
            %         L_b_initial(:,i)=L_b_initial(:,i)/sum(L_b_initial(:,i));
            
            
            
            %%%%%%%%%%
%             tempo=zeros(2^r,1);
%             tempo(randi(2^r))=1;
%             tempo=rand(2^r,1)/20+tempo;
%             tempo=tempo./sum(tempo);
%             L_b(:,i)=tempo;
            
            
        end
    end
    
    for i=1:1:n
        R_b(:,i)=vector_to_cat(R_i(:,i)',r);
        
        if rand<perturbation_prob
            %%% slight perturnation
            R_b(:,i)=R_b(:,i)+1/2^r;
            R_b(:,i)=R_b(:,i)/sum(R_b(:,i));
            
            %%%% replacing with random vectors
            %         R_b_initial(:,i)=rand(2^r,1);
            %         R_b_initial(:,i)=R_b_initial(:,i)/sum(R_b_initial(:,i));
            
            %%%%%%
%             tempo=zeros(2^r,1);
%             tempo(randi(2^r))=1;
%             tempo=rand(2^r,1)/20+tempo;
%             tempo=tempo./sum(tempo);
%             R_b(:,i)=tempo;
        end
        
    end
    
    
    %%%%%%%%%  uniform initialization
    
    %     L_b=ones(2^r,m)/2^r;
    %
    %     R_b=ones(2^r,n)/2^r;
    
    
    
    
    [err, A_recovered, L_b, R_b]=MyAlgorithm_feeded( original(:,:,sample),masks(:,:,sample),m,n,r,"second", L_b, R_b, C0, C1);
    
    err
    my_errors(sample)=min(my_errors(sample),err);
 
    
end
end

end












function [error, A_recovered, L_b, R_b]=MyAlgorithm_feeded(A,obs_pattern,m,n,r,mode, Lb_init, Rb_init, C0, C1)
%function [error]=MyAlgorithm()



% m=100;
% n=100;
% r=5;
%p_obs=.35;
iterations=50;
error_threshold=10^-6; % for BP convergence
% m,n=size(A);
% L_i=randi(2,m,r)-1;
% R_i=randi(2,r,n)-1;
% A=mod(L_i*R_i,2);
%
%
%
% %%%%Random observation pattern
% obs_pattern=rand(m,n);
% obs_pattern=(obs_pattern<p_obs);



L_connections=cell(1,m);
L_degree_profile=cell(1,m);
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
R_degree_profile=cell(1,n);
for j=1:n
    temp=[];
    for i=1:m
        if obs_pattern(i,j)==1
            temp=[temp,i];
        end
    end
    R_connections{j}=temp;
end






%%% Initializing near the ground truth

%
% L_b_initial=zeros(2^r,m)/2^r;
% R_b_initial=zeros(2^r,n)/2^r;
% perturbation_prob=.95;
% for i=1:1:m
%     L_b_initial(:,i)=vector_to_cat(L_i(i,:),r);
%
%     %%% perturbing from the solution
%     if rand<perturbation_prob
% %         L_b_initial(:,i)=L_b_initial(:,i)+rand(2^r,1)/20;
%         L_b_initial(:,i)=rand(2^r,1);
%         L_b_initial(:,i)=L_b_initial(:,i)/sum(L_b_initial(:,i));
%     end
%
%
% end
%
% for i=1:1:n
%     R_b_initial(:,i)=vector_to_cat(R_i(:,i)',r);
%
%     %%% perturbing from the solution
%     if rand<perturbation_prob
%         %R_b_initial(:,i)=R_b_initial(:,i)+rand(2^r,1)/20;
%         R_b_initial(:,i)=rand(2^r,1);
%         R_b_initial(:,i)=R_b_initial(:,i)/sum(R_b_initial(:,i));
%     end
% end




if mode=="first"
    
    L_b_initial=ones(2^r,m)/2^r;
    
    R_b_initial=ones(2^r,n)/2^r;
    
end

if mode=="second"
    
    L_b_initial=Lb_init;
    
    %     L_b_initial=Lb_init+ones(2^r,m)/2^r;
    %     L_b_initial(:,i)=L_b_initial(:,i)/sum(L_b_initial(:,i));
    
    R_b_initial=Rb_init;
    %
    %
    %     R_b_initial=Rb_init+ones(2^r,n)/2^r;
    %     R_b_initial(:,i)=R_b_initial(:,i)/sum(R_b_initial(:,i));
    
end






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
    i;
    %     L_b_initial=rand(2^r,m);
    %     L_b_initial=L_b_initial./sum(L_b_initial)
    %
    %     R_b_initial=rand(2^r,n);
    %     R_b_initial=R_b_initial./sum(R_b_initial);
    
    [converged,L_b,R_b]= Recovery(m,n,r,iterations,error_threshold,A, obs_pattern,C0,C1,L_b_initial,R_b_initial,L_b_fixed_vals);
    convergence_status=[convergence_status,converged];%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    L_b(:,i);
    L_b;
    R_b;
    fixed_value=discretesample(L_b(:,i), 1);
    L_b_fixed_vals=[L_b_fixed_vals,double(1:2^r == fixed_value)'];
    L_b_fixed_vals;
    
    
    nnz(L_b);%%%%%%%%%%%%%%%%%%%
    nnz(R_b);%%%%%%%%%%%%%%%%%
    if nnz(L_b)<=m && nnz(R_b)<=n
        sum(L_b);
        sum(R_b);
        break
    end
    
    
    
end



L=zeros(m,r);

for i=1:1:m
    loc=find(L_b(:,i));
    
    if isempty(loc)
        loc=randi(2^r);
    end
    
    L(i,:)=de2bi(loc(1)-1,r);
end

R=zeros(r,n);


for i=1:1:n
    loc=find(R_b(:,i));
    if isempty(loc)
        loc=randi(2^r);
    end
    R(:,i)=de2bi(loc(1)-1,r)';
end


convergence_status;%%%%%%%%%%%%%%%%



A_recovered=(L*R)>0;
%error=nnz((A-A_recovered).*obs_pattern)/nnz(obs_pattern);
error=nnz(A-A_recovered)/(m*n);


end

function[out]=vector_to_cat(a,r)

id=bi2de(a);
out=zeros(2^r,1);
out(id+1)=1;

end

function[converged,L_b,R_b]= Recovery(m,n,r,iterations,error_threshold,A, obs_pattern,C0,C1,L_b,R_b,L_b_fixed_vals)




%%%%%%%%%% everytime starts from uniform beleif

% L_b=ones(2^r,m)/2^r;
% R_b=ones(2^r,n)/2^r;




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
                L_b_new(:,i)=temp/sum(temp).*L_b_new(:,i);
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




















