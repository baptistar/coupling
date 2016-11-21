function [model_mean, model_cov, model_out] = gaussian_linearization(particle_system, model)
% GAUSSIAN_LINEARIZATION: Function computes the mean and covariance of the
% Gaussian linearization to the non-linear model by running the model at the
% linearization point.
%
% See COPYRIGHT for license information 
%

% Extract inputs from particle_system
qoi             = particle_system.qoi;
eng_problem     = particle_system.eng_problem;

% Run engineering problem at the given model
[dR_dx, dR_dy, input_mean, input_cov, model_out, output_vars] = eng_problem(particle_system, model);

% Declare function for searching cell contents
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));

%% Find and Extract QoI
n_qoi = length(qoi);
qoi_idx = [];
for i=1:n_qoi
    qoi_idx = [qoi_idx, find(cellfun(cellfind(qoi{i}),output_vars))];
end

%% Find Mean

% Setup mean with model output for specified QoI
mean_orig  = model_out.all_vars;
model_mean = mean_orig(qoi_idx);

%% Find Covariance

% Assemble covariance matrix and extract specified QoI
model_cov = (dR_dy\dR_dx)*input_cov*(dR_dy\dR_dx)';
model_cov = model_cov(qoi_idx, qoi_idx);

end
