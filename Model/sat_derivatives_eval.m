function [dR_dx_eval, dR_dy_eval] = sat_derivatives_eval(var_info, model)
% SAT_DERIVATIVES: Function computes the derivatives for the satellite
% model based on symbolic differentiation at the given inputs in
% variable_info.
%
% See COPYRIGHT for license information 
%
% The satellite problem is adapted from the reference:
% Sankararaman, S., and Mahadevan, S. Likelihood-based approach to multidisciplinary
% analysis under uncertainty. Journal of Mechanical Design, 2012.

% Extract variables
if exist('Model/sat_derivatives_data.mat','file')
	load('Model/sat_derivatives_data','dR_dx','dR_dy');
else
	[dR_dx, dR_dy] = sat_derivatives_setup(var_info, length(model));
end

%% Setup Values

% Determine the number of inputs and outputs
n_rand_inputs = length(var_info.rand_inputs);
n_outputs     = length(var_info.output_vars);
n_residuals   = size(dR_dy,1);

% Set values for fixed inputs variables
RE         = var_info.RE;
mu         = var_info.mu;
phi        = var_info.phi;
c          = var_info.c;
A_s        = var_info.A_s;
sun_i      = var_info.sun_i;
dt_slew    = var_info.dt_slew;
M          = var_info.M;
rho        = var_info.rho;
A          = var_info.A;
n          = var_info.n;
omega_max  = var_info.omega_max;
P_hold     = var_info.P_hold;
I_d        = var_info.I_d;
nu         = var_info.nu;
LT         = var_info.LT;
eps_deg    = var_info.eps_deg;
r_lw       = var_info.r_lw;
n_sa       = var_info.n_sa;
rho_sa     = var_info.rho_sa;
t          = var_info.t;
D          = var_info.D;
I_bodyx    = var_info.I_bodyx;
I_bodyy    = var_info.I_bodyy;
I_bodyz    = var_info.I_bodyz;

% Set values for random inputs variables
H          = var_info.H;
Po         = var_info.Po;
F_s        = var_info.F_s;
theta      = var_info.theta;
L_sp       = var_info.L_sp;
q          = var_info.q;
RD         = var_info.RD;
L_a        = var_info.L_a;
C_d        = var_info.C_d;

% Set values for orbit output variables
v          = var_info.v;
dt_orbit   = var_info.dt_orbit;
dt_eclipse = var_info.dt_eclipse;
theta_slew = var_info.theta_slew;

% Set values for attitude control output variables
tau_slew   = var_info.tau_slew;
tau_g      = var_info.tau_g;
tau_sp     = var_info.tau_sp;
tau_m      = var_info.tau_m;
tau_a      = var_info.tau_a;
tau_dist   = var_info.tau_dist;
tau_tot    = var_info.tau_tot;
PACS       = var_info.PACS;
P_tot      = var_info.P_tot;

% Set values for power output variables
T_e        = var_info.T_e;
T_d        = var_info.T_d;
P_e        = var_info.P_e;
P_d        = var_info.P_d;
P_sa       = var_info.P_sa;
PBOL       = var_info.PBOL;
PEOL       = var_info.PEOL;
A_sa       = var_info.A_sa;

% Compute solar array geometry quantities
L          = var_info.L;
W          = var_info.W;
m_sa       = var_info.m_sa;

% Compute solar array moments of inertia
I_sax      = var_info.I_sax;
I_say      = var_info.I_say;
I_saz      = var_info.I_saz;
I_max      = var_info.I_max;
I_min      = var_info.I_min;

%% Evaluate Derivatives Numerically

% Generate matrices to store derivatives
dR_dx_eval = zeros(n_residuals,n_rand_inputs);
dR_dy_eval = zeros(n_residuals,n_outputs);

% Evaluate derivatives with respect to inputs
for i=1:n_residuals
    for j=1:n_rand_inputs
        dR_dx_eval(i,j) = eval(subs(dR_dx{i,j})); 
    end
end

% Evaluate derivatives with respect to outputs
for i=1:n_residuals
    for j=1:n_outputs
        dR_dy_eval(i,j) = eval(subs(dR_dy{i,j})); 
    end
end

% -- END OF FILE --