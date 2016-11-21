function [dR_dx, dR_dy] = sat_derivatives_setup(var_info, n_vars)
% SAT_DERIVATIVES_SETUP: Function setups the derivative functions for the
% satellite model based on symbolic differentiation at the given inputs in
% var_info.
%
% See COPYRIGHT for license information 
%
% The satellite problem is adapted from the reference:
% Sankararaman, S., and Mahadevan, S. Likelihood-based approach to multidisciplinary
% analysis under uncertainty. Journal of Mechanical Design, 2012.

%% Run model

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

% Compute outputs
var_info = run_sat_model(var_info, ones(1,n_vars));

%% Define Variables

% Define all variables
rand_inputs   = var_info.rand_inputs;
fixed_inputs  = var_info.fixed_inputs;
output_vars   = var_info.output_vars;

% Determine the number of inputs and outputs
n_rand_inputs = length(rand_inputs);
n_fix_inputs  = length(fixed_inputs);
n_outputs     = length(output_vars);

% Convert variables to symbolic form
for i=1:n_rand_inputs
    syms(rand_inputs{i});
end

for i=1:n_fix_inputs
    syms(fixed_inputs{i});
end

for i=1:n_outputs
    syms(output_vars{i});
end

%% Define Residual Functions

% Define residual functions for orbit subsystem
R_1 = @(v, RE, H, mu) v - sqrt(mu/(RE + H));
R_2 = @(dt_orbit, RE, H, v) dt_orbit - 2*pi*(RE + H)/v;
R_3 = @(dt_eclipse, RE, H, dt_orbit) dt_eclipse - dt_orbit/pi*asin(RE/(RE+H));
R_4 = @(theta_slew, RE, H, phi) theta_slew - atan(sin(phi/RE)/(1 - cos(phi/RE) + H/RE));

% Define residual functions for attitude control subsystem
R_5  = @(tau_slew, theta_slew, dt_slew, I_max) tau_slew - 4*theta_slew/(dt_slew)^2*I_max;
R_6  = @(tau_g, RE, H, mu, v, theta, I_max, I_min) tau_g - 3*mu/(2*(RE+H)^3)*abs(I_max - I_min)*sin(2*theta);
R_7  = @(tau_sp, L_sp, F_s, c, A_s, q, sun_i) tau_sp - L_sp*F_s/c*A_s*(1+q)*cos(sun_i);
R_8  = @(tau_m, M, RD, RE, H) tau_m - 2*M*RD/(RE+H)^3;
R_9  = @(tau_a, L_a, rho, C_d, A, v) tau_a - 0.5*L_a*rho*C_d*A*v^2;
R_10 = @(tau_dist, tau_g, tau_sp, tau_m, tau_a) tau_dist - sqrt(tau_g^2 + tau_sp^2 + tau_m^2 + tau_a^2);
if var_info.tau_slew >= var_info.tau_dist
    R_11 = @(tau_tot, tau_slew) tau_tot - tau_slew;
elseif var_info.tau_dist > var_info.tau_slew
    R_11 = @(tau_tot, tau_dist) tau_tot - tau_dist;
end
R_12 = @(PACS, tau_tot, omega_max, n, P_hold) PACS - (tau_tot*omega_max + n*P_hold);

% Define residual functions for power subsystem
R_13 = @(P_tot, PACS, Po) P_tot - (PACS + Po);
R_14 = @(T_e, dt_eclipse) T_e - dt_eclipse;
R_15 = @(T_d, dt_orbit, T_e) T_d - (dt_orbit - T_e);
R_16 = @(P_e, P_tot) P_e - P_tot;
R_17 = @(P_d, P_tot) P_d - P_tot;
R_18 = @(P_sa, P_e, T_e, P_d, T_d) P_sa - (P_e*T_e/0.6 + P_d*T_d/0.8)/T_d;
R_19 = @(PBOL, nu, F_s, I_d, sun_i) PBOL - nu*F_s*I_d*cos(sun_i);
R_20 = @(PEOL, PBOL, eps_deg, LT) PEOL - PBOL*(1-eps_deg)^(LT);
R_21 = @(A_sa, P_sa, PEOL) A_sa - P_sa/PEOL;

% Compute solar array geometry quantities
R_22 = @(L, A_sa, r_lw, n_sa) L - sqrt(A_sa*r_lw/n_sa);
R_23 = @(W, A_sa, r_lw, n_sa) W - sqrt(A_sa/(r_lw*n_sa));
R_24 = @(m_sa, rho_sa, L, W, t) m_sa - 2*rho_sa*L*W*t;

% Compute solar array moments of inertia
R_25 = @(I_sax, m_sa, L, D, t) I_sax - m_sa*(1/12*(L^2 + t^2) + (D + L/2)^2);
R_26 = @(I_say, m_sa, t, W) I_say - m_sa/12*(t^2 + W^2);
R_27 = @(I_saz, m_sa, L, D, W) I_saz - m_sa*(1/12*(L^2 + W^2) + (D + L/2)^2);
if ((var_info.I_sax + var_info.I_bodyx) == max([var_info.I_sax, var_info.I_say, var_info.I_saz] + [var_info.I_bodyx, var_info.I_bodyy, var_info.I_bodyz]))
    R_28 = @(I_max, I_sax, I_bodyx) I_max - (I_sax + I_bodyx);
elseif ((var_info.I_say + var_info.I_bodyy) == max([var_info.I_sax, var_info.I_say, var_info.I_saz] + [var_info.I_bodyx, var_info.I_bodyy, var_info.I_bodyz]))
    R_28 = @(I_max, I_say, I_bodyy) I_max - (I_say + I_bodyy);
elseif ((var_info.I_saz + var_info.I_bodyz) == max([var_info.I_sax, var_info.I_say, var_info.I_saz] + [var_info.I_bodyx, var_info.I_bodyy, var_info.I_bodyz]))
    R_28 = @(I_max, I_saz, I_bodyz) I_max - (I_saz + I_bodyz);
end
if ((var_info.I_sax + var_info.I_bodyx) == min([var_info.I_sax, var_info.I_say, var_info.I_saz] + [var_info.I_bodyx, var_info.I_bodyy, var_info.I_bodyz]))
    R_29 = @(I_min, I_sax, I_bodyx) I_min - (I_sax + I_bodyx);
elseif ((var_info.I_say + var_info.I_bodyy) == min([var_info.I_sax, var_info.I_say, var_info.I_saz] + [var_info.I_bodyx, var_info.I_bodyy, var_info.I_bodyz]))
    R_29 = @(I_min, I_say, I_bodyy) I_min - (I_say + I_bodyy);
elseif ((var_info.I_saz + var_info.I_bodyz) == min([var_info.I_sax, var_info.I_say, var_info.I_saz] + [var_info.I_bodyx, var_info.I_bodyy, var_info.I_bodyz]))
    R_29 = @(I_min, I_saz, I_bodyz) I_min - (I_saz + I_bodyz);
end

% Assemble residual_funcs
residual_funcs = {R_1, R_2, R_3, R_4, R_5, R_6, R_7, R_8, R_9, R_10, ...
                  R_11, R_12, R_13, R_14, R_15, R_16, R_17, R_18, R_19, R_20, ...
                  R_21, R_22, R_23, R_24, R_25, R_26, R_27, R_28, R_29};

% Compute number of residual functions
n_residuals = length(residual_funcs);

%% Evaluate Derivatives Symbolically

% Generate cells to store derivatives
dR_dx = cell(n_residuals,n_rand_inputs);
dR_dy = cell(n_residuals,n_outputs);

% Compute derivatives with respect to random inputs symbolically
for i=1:n_residuals
    for j=1:n_rand_inputs
        dR_dx{i,j} = diff(residual_funcs{i}, sym(rand_inputs{j})); 
    end
end

% Compute derivatives with respect to outputs symbolically
for i=1:n_residuals
    for j=1:n_outputs
        dR_dy{i,j} = diff(residual_funcs{i}, sym(output_vars{j})); 
    end
end

% Save dR_dx and dR_dy
save('Model/sat_derivatives_data','dR_dx','dR_dy');

% -- END OF FILE --