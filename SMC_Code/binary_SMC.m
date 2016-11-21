function opt_model = binary_SMC(input)
% BINARY_SMC: Function explores the binary space of models and returns the
% expected value for each binary variable.
%
% See COPYRIGHT for license information 
%

% Declare recommended algorithm parameters
init_rho     = 1e-5;
alpha_star   = 0.90;

%% Setup Parameters

% Setup struct
ps = struct;

% Save parameters in struct
ps.n_models    = input.n_models;
ps.lambda      = input.lambda;
ps.qoi         = input.qoi;
ps.n_vars      = input.n_vars;
ps.eng_problem = input.eng_problem;

% Find the reference model mean & covariance
[ps.ref_mean, ps.ref_cov, ps.ref_mean_out] = gaussian_linearization(ps, ones(1,ps.n_vars));

%% Setup SMC sampler

% Save initial eps
ps.rho = init_rho;

% Save alpha
ps.alpha_star = alpha_star;

% Initialize mh_accept
ps.mh_accept = 1;

% Sample initial models
ps.models = sample_models(ps.n_models, ps.n_vars);

% Compute model distances and model_prob
ps.model_val = model_objective(ps, ps.models);

% Set initial importance weights
ps.weights = 1/ps.n_models*ones(ps.n_models,1);

% Initialize logistic_A
ps.logistic_A = eye(ps.n_vars);

% Initialize counter
counter = 1;

% Setup cells to store models and weights
models_iter  = {ps.models};
weights_iter = {ps.weights};

% Set initial alpha and weights
ps = find_tolerance(ps);

% Setup vectors to store ESS and PDIV
ps.ess  = eff_sample_size(ps.weights);
ps.pdiv = unique_particles(ps.models);

%% Run SMC sampler

disp('-----------------------------------------------------');
disp('Binary Sampling with SMC');
disp('-----------------------------------------------------');

% Print Current Results
fprintf('Iteration %3d, EPS %5.3e, ESS %5.1f, UNIQUE %4.3f, MH ACCEPT: %4.3f\n', ...
        counter, ps.rho, ps.ess, ps.pdiv, ps.mh_accept);

% Run algorithm until convergence of rho_k
while(ps.pdiv > input.delta)

    % Fit logistic regression model
    ps = fit_logistic(ps);

	% Resample models
    ps = smc_resample(ps);

    % Move particles with the logistic regression model
    ps = move_particles(ps);
    
    % Find step length and compute weights
    ps = find_tolerance(ps);
    
    % Update counter
    counter = counter + 1;

    % Update ESS and PDIV
    ps.ess  = eff_sample_size(ps.weights);
    ps.pdiv = unique_particles(ps.models);

    % Add current models to cell
    models_iter{counter}  = ps.models;
    weights_iter{counter} = ps.weights;

    % Print Current Results
    fprintf('Iteration %3d, EPS %5.3e, ESS %5.1f, UNIQUE %4.3f, MH ACCEPT: %4.3f\n', ...
            counter, ps.rho, ps.ess, ps.pdiv, ps.mh_accept);
                
end

% Extract models from final iteration
final_models    = ps.models;
final_model_val = ps.model_val;

% Find unique models
[unique_models, unique_idx] = unique(final_models,'rows');
unique_model_val = final_model_val(unique_idx);

% Determine optimal model
opt_idx = find(unique_model_val == max(unique_model_val));
opt_model = unique_models(opt_idx,:);

% Save Results
lambda_str = strrep(num2str(input.lambda,'%.1e'),'.','p');
save(['Results/smc_results_' lambda_str '_lambda.mat'])

% -- END OF FILE --
