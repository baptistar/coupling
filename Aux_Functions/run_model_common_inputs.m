function [mc_input, output_cell] = run_model_common_inputs(ps, models, n_samples)
% RUN_MODEL_COMMON_INPUTS: Function runs Monte Carlo simulations on the original model 
%
% See COPYRIGHT for license information 
%

%% Determine Sample inputs
mc_input  = zeros(n_samples, length(ps.var_info.rand_inputs));

for i=1:n_samples

    % Setup random variables
    H       = ps.var_info.H_std*randn + ps.var_info.H_mean;
    Po      = ps.var_info.Po_std*randn + ps.var_info.Po_mean;
    F_s     = ps.var_info.F_s_std*randn + ps.var_info.F_s_mean;
    theta   = ps.var_info.theta_std*randn + ps.var_info.theta_mean;
    L_sp    = ps.var_info.L_sp_std*randn + ps.var_info.L_sp_mean;
    q       = ps.var_info.q_std*randn + ps.var_info.q_mean;
    RD      = ps.var_info.RD_std*randn + ps.var_info.RD_mean;
    L_a     = ps.var_info.L_a_std*randn + ps.var_info.L_a_mean;
    C_d     = ps.var_info.C_d_std*randn + ps.var_info.C_d_mean;

    % Save inputs
    mc_input(i,:) = [H, Po, F_s, theta, L_sp, q, RD, L_a, C_d];

end

% Declare function for searching cell contents
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));

% Determine QoI indices
n_qoi = length(ps.qoi);
qoi_idx = [];
for i=1:n_qoi
    qoi_idx = [qoi_idx, find(cellfun(cellfind(ps.qoi{i}), ps.var_info.output_vars))];
end

% Declare cell to store output from each model
n_models = size(models,1);
output_cell = cell(n_models,1);

for j=1:n_models

    fprintf('Running MC on Model %d/%d\n', j, n_models);

    % Declare matrix to store inputs and outputs
    mc_output = zeros(n_samples, length(ps.var_info.output_vars));

    % Compute and save outputs for each input
    for i=1:n_samples

        ps.var_info.H       = mc_input(i,1);
        ps.var_info.Po      = mc_input(i,2);
        ps.var_info.F_s     = mc_input(i,3);
        ps.var_info.theta   = mc_input(i,4);
        ps.var_info.L_sp    = mc_input(i,5);
        ps.var_info.q       = mc_input(i,6);
        ps.var_info.RD      = mc_input(i,7);
        ps.var_info.L_a     = mc_input(i,8);
        ps.var_info.C_d     = mc_input(i,9);

        options = optimoptions('fsolve','Display','none','TolFun',1e-15);
        mc_output(i,:) = fsolve(@(out) satellite_solver(out, ps.var_info, models(j,:)), rand(29,1), options);

    end

    % Extract components of QoI from mc_output and save in output_cell
    output_cell{j} = mc_output(:,qoi_idx);

end

end
