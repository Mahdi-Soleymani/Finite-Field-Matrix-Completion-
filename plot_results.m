function [out]=plot_results()

num_of_iterations=30;
 p_obs=.15:.05:.25;
%p_obs=.1:.05:.4;
%p_obs=1;
accuracy_bit=zeros(length(p_obs),1);
accuracy_block=zeros(length(p_obs),1);

out=zeros(length(p_obs),1);

for k=1:length(p_obs)
    
for i=1:num_of_iterations
    p_obs(k)
    i
    [error]=MyAlgorithm_boolean_multiplication_max_probbability_assignement(p_obs(k));
    %[error]=MyAlgorithm_boolean_multiplication_stopping_condiition(p_obs(k))
    %[error]=MyAlgorithm_boolean_multiplication_picked_by_entropy(p_obs(k));
    %[error]=MyAlgorithm_boolean_multiplication(p_obs(k));
% [error]=MyAlgorithm(p_obs(k));
%[error]=MyAlgorithm_fix_higher_degrees_first(p_obs(k));


if error==0
    accuracy_block(k)=accuracy_block(k)+1;
end

accuracy_bit(k)=accuracy_bit(k)+error;

end
end
out=1-accuracy_bit/num_of_iterations;
accuracy_bit
accuracy_block
plot(p_obs,out);


end