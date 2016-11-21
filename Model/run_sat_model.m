function [var_info] = run_sat_model(var_info, model)
% RUN_SAT_MODEL: Function runs satellite model for the given coupling
% variables and saves new outputs in var_info.
%
% See COPYRIGHT for license information 
%
% The satellite problem is adapted from the reference:
% Sankararaman, S., and Mahadevan, S. Likelihood-based approach to multidisciplinary
% analysis under uncertainty. Journal of Mechanical Design, 2012.

% Compute outputs
options = optimoptions('fsolve','Display','none','TolFun',1e-15);
outputs = fsolve(@(out) satellite_solver(out, var_info, model), rand(29,1), options);

% Save output variables
var_info.v          = outputs(1);
var_info.dt_orbit   = outputs(2);
var_info.dt_eclipse = outputs(3);
var_info.theta_slew = outputs(4);
var_info.tau_slew   = outputs(5);
var_info.tau_g      = outputs(6);
var_info.tau_sp     = outputs(7);
var_info.tau_m      = outputs(8);
var_info.tau_a      = outputs(9);
var_info.tau_dist   = outputs(10);
var_info.tau_tot    = outputs(11);
var_info.PACS       = outputs(12);
var_info.P_tot      = outputs(13);
var_info.T_e        = outputs(14);
var_info.T_d        = outputs(15);
var_info.P_e        = outputs(16);
var_info.P_d        = outputs(17);
var_info.P_sa       = outputs(18);
var_info.PBOL       = outputs(19);
var_info.PEOL       = outputs(20);
var_info.A_sa       = outputs(21);
var_info.L          = outputs(22);
var_info.W          = outputs(23);
var_info.m_sa       = outputs(24);
var_info.I_sax      = outputs(25);
var_info.I_say      = outputs(26);
var_info.I_saz      = outputs(27);
var_info.I_max      = outputs(28);
var_info.I_min      = outputs(29);

% Save all values
var_info.all_vars = outputs;

end
