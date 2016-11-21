function model_val = model_objective(particle_system, model)
% MODEL_OBJECTIVE: Function evaluates the objective function for the
% specified binary model based on the applied prior and KL divergence.
%
% See COPYRIGHT for license information 
%

% Extract variables from particle_system
n_vars     = particle_system.n_vars;
lambda     = particle_system.lambda;

% Determine the number of models
[n_models, ~] = size(model);

% Declare vector to store model_val
model_val = zeros(n_models,1);

% Setup un-normalized prior for the number of binary connections in the model
coupling_prob    = n_vars+1:-1:1;
coupling_prior   = coupling_prob;

for i=1:n_models

    % Evaluate Prior
    model_prior = lambda*coupling_prior(sum(model(i,:)) + 1);

    % Evaluate KL Divergence based on Gaussian linearization
    model_dist = kldiv_likelihood(particle_system, model(i,:));
    
    % Evaluate negative of objective function (maximize objective)
    model_val(i) = -(model_dist - model_prior);

end

end
