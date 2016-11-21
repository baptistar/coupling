function [binary_models] = sample_models(n_models, n_vars)
% SAMPLE_MODELS: Function samples the binary models to use in the SMC
% algorithm based on an l1-norm prior penalty on the model space.
%
% See COPYRIGHT for license information 
%

% Generate matrix of zeros with ones along diagonals
binary_models = zeros(n_models,n_vars);

% Setup prior for the number of binary connections in the model
coupling_prob    = n_vars+1:-1:1;
coupling_nmodels = pascal_row(n_vars);
coupling_prior   = coupling_prob/sum(coupling_prob.*coupling_nmodels);

% Sample each binary vector individually
for i=1:n_models
    
    % Sample the number of elements in the model
    n_coupling = discrete_sample(coupling_prior,1)-1;

    % Sample mask entries based on n_coupling
    loc_entry  = datasample(1:n_vars, n_coupling, 'Replace', false);
    
    % Save entries in binary_models
    binary_models(i,loc_entry) = 1;
    
end

end