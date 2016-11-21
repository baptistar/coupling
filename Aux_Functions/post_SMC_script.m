% POST SMC Script
%
% See COPYRIGHT for license information 
%

clear; close all; clc

addpath SMC_Code
addpath Model
addpath Aux_Functions
addpath Results

% Setup lambda vector
lambda_folders = dir('Results/*_Lambda');
n_folders = length(lambda_folders);

% Setup parameters
num_mc_samples   = 10000;

%% Compute Gaussian approximation and MC samples

% Setup cell/vector to store optimal model information
lambda_vect         = zeros(n_folders,1);
opt_models_cell     = cell(n_folders,1);
opt_models_kl_div   = zeros(n_folders,1);
opt_models_coupls   = zeros(n_folders,1);
opt_model_gaussian  = cell(n_folders,2);

for i=1:n_folders

    fprintf('Analyzing Model: %d/%d\n', i, n_folders);

    % Load results in data file
    load(['Results/' lambda_folders(i).name '/smc_results'])
    
    % Extract models from final iteration
    final_models    = ps.models;
    final_model_val = ps.model_val;

    % Find unique models
    [unique_models, unique_idx] = unique(final_models,'rows');
    unique_model_val = final_model_val(unique_idx);

    % Determine optimal model (or optimal models if multiple models have the 
    % same value for the objective function)
    opt_idx = find(unique_model_val == max(unique_model_val));
    opt_models_cell{i} = unique_models(opt_idx,:);

    % Compute KL Divergence and sparsity for these models
    opt_models_kl_div(i) = kldiv_likelihood(ps, opt_models_cell{i});
    opt_models_coupls(i) = sum(opt_models_cell{i});

    % Run Gaussian approximation on model
    [dec_mean, dec_cov, ~] = gaussian_linearization(ps, opt_models_cell{i});
    opt_model_gaussian{i,1}  = dec_mean; 
    opt_model_gaussian{i,2}  = dec_cov;

    % Save lambda value
    lambda_vect(i) = ps.lambda;

end

% Find unique set of models
[all_models, unique_idx] = unique(cell2mat(opt_models_cell),'rows');

% Find unique set of values for post-processed outputs
lambda_str         = {lambda_folders.name};
lambda_str         = lambda_str(unique_idx);
lambda_vect        = lambda_vect(unique_idx);
opt_models_kl_div  = opt_models_kl_div(unique_idx);
opt_models_coupls  = opt_models_coupls(unique_idx);
opt_model_gaussian = opt_model_gaussian(unique_idx,:);

% Append reference to all models
ref_model   = ones(1,ps.n_vars);
all_models  = [all_models; ref_model];

% Find unique set of models
[all_models, unique_idx] = unique(all_models,'rows');

% Run full Monte Carlo analysis on all models
[model_mc_full_input, model_mc_full_output] = run_model_common_inputs(ps, all_models, num_mc_samples);

ref_model_mc_full_output = model_mc_full_output{end};
model_mc_full_output     = model_mc_full_output(1:end-1);

% Save data
save(['Results/post_SMC_data'])

% -- END OF FILE --
