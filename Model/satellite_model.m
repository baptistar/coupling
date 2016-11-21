function [dR_dx, dR_dy, input_mean, input_cov, model_out, output_vars] = satellite_model(particle_system, model)
% SATELLITE_MODEl_DERS: Function evaluates the derivatives of the satellite model and the 
% first-order mean of the system. The function also returns the mean and covariance of the input
% variables and the output_vars in order to compare to the specified QoI
%
% See COPYRIGHT for license information 
%
% The satellite problem is adapted from the reference:
% Sankararaman, S., and Mahadevan, S. Likelihood-based approach to multidisciplinary
% analysis under uncertainty. Journal of Mechanical Design, 2012.

% Setup problem parameters

% Declare disciplines, coupling variables and variable data
disciplines     = {'orbit','ac','power'};
model_couplings = {'power.dt_orbit','power.dt_eclipse','power.PACS', ...
				   'ac.I_max','ac.I_min','ac.dt_orbit','ac.v','ac.theta_slew'};
decoupling_idx  = [1,5,13,30];

% Setup var_info
var_info = satellite_setup(particle_system);

% Declare values for random variables
var_info.H          = var_info.H_mean;
var_info.Po         = var_info.Po_mean;
var_info.F_s        = var_info.F_s_mean;
var_info.theta      = var_info.theta_mean;
var_info.L_sp       = var_info.L_sp_mean;
var_info.q          = var_info.q_mean;
var_info.RD         = var_info.RD_mean;
var_info.L_a        = var_info.L_a_mean;
var_info.C_d        = var_info.C_d_mean;

%% Compute Output Values and Derivatives

% Compute outputs
var_info = run_sat_model(var_info, model);

% Compute derivatives
[dR_dx, dR_dy] = sat_derivatives_eval(var_info, model);

%% Decouple Disciplines (Setup Mask for dR_dy)

% Declare function for searching cell contents
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));

% Find decouple_vars indices
decouple_vars = find(~model);

% Set indices corresponding to decoupled variables in dR_dy_eval to 0
for k=1:length(decouple_vars)
    
    % Extract coupling to mask
    coupling_elim = model_couplings{decouple_vars(k)};
    
    % Seperate variable and discipline
    strings = strsplit(coupling_elim,'.');
    coupling_dis = strings{1};
    coupling_var = strings{2};
    
    % Find discipline and coupling variable indices
    discp_index = find(cellfun(cellfind(coupling_dis),disciplines));
    var_index = find(cellfun(cellfind(coupling_var),var_info.output_vars));
    
    % Find discipline_indices
    mask_index = decoupling_idx(discp_index):(decoupling_idx(discp_index+1)-1);
    
    % Mask dR_dy
    dR_dy(mask_index, var_index) = 0;
    
end

%% Setup Input Parameters

% Compute input_mean
input_mean = [var_info.H_mean, var_info.Po_mean, var_info.F_s_mean, ...
              var_info.theta_mean, var_info.L_sp_mean, var_info.q_mean, ...
              var_info.RD_mean, var_info.L_a_mean, var_info.C_d_mean];

% Compute input_cov
input_cov = diag([var_info.H_std^2, var_info.Po_std^2, var_info.F_s_std^2, ...
				  var_info.theta_std^2, var_info.L_sp_std^2, var_info.q_std^2, ...
				  var_info.RD_std^2, var_info.L_a_std^2, var_info.C_d_std^2]);

% Save model out & output_vars
model_out   = var_info;
output_vars = var_info.output_vars;

end