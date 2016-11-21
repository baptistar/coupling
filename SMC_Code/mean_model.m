function mean_vector = mean_model(particle_system)
% MEAN_MODEL: Function computes mean binary model based on model weights.
%
% See COPYRIGHT for license information 
%

% Extract model inputs
models  = particle_system.models;
weights = particle_system.weights;

% Compute n_vars
n_vars = size(models,2);

% Declare model mean vector
mean_vector = zeros(n_vars,1);

% Compute mean vector components
for i=1:n_vars
    mean_vector(i) = models(:,i)'*weights;
end