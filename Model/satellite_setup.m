function var_info = satellite_setup(particle_system)
% SATELLITE_SETUP: Function returns a struct containing the values for the
% fixed inputs (deterministic values) and random input parameters (mean 
% and covariance).
%
% See COPYRIGHT for license information 
%
% The satellite problem is adapted from the reference:
% Sankararaman, S., and Mahadevan, S. Likelihood-based approach to multidisciplinary
% analysis under uncertainty. Journal of Mechanical Design, 2012.

%% Declare parameters

% Define variables names
rand_inputs  = {'H','Po','F_s','theta','L_sp','q','RD','L_a','C_d'};
fixed_inputs = {'RE','mu','phi','c','A_s','sun_i','dt_slew','M','rho', ...
                'A','n','omega_max','P_hold','I_d','nu','LT','eps_deg',...
                'r_lw','n_sa','rho_sa','t','D','I_bodyx','I_bodyy','I_bodyz'};
output_vars  = {'v','dt_orbit','dt_eclipse','theta_slew','tau_slew', ...
                'tau_g', 'tau_sp', 'tau_m', 'tau_a', 'tau_dist', 'tau_tot', ...
                'PACS', 'P_tot', 'T_e', 'T_d', 'P_e', 'P_d',  ...
                'P_sa', 'PBOL', 'PEOL', 'A_sa', 'L', 'W', 'm_sa',  ...
                'I_sax', 'I_say', 'I_saz', 'I_max', 'I_min'};


% Declare fixed quantities
RE         = 6378140;
mu         = 3.986*10^14;
phi        = 235;
c          = 2.9989*10^8;
A_s        = 13.85;
sun_i      = 0;
dt_slew    = 760;
M          = 7.96*10^15;
rho        = 5.1480*10^(-11);
A          = 13.85;
n          = 3;
omega_max  = 6000;
P_hold     = 20;
I_d        = 0.77;
nu         = 0.22;
LT         = 15;
eps_deg    = 0.0375;
r_lw       = 3;
n_sa       = 3;
rho_sa     = 700;
t          = 0.005;
D          = 2;
I_bodyx    = 6200;
I_bodyy    = 6200;
I_bodyz    = 4700;

% Declare random quantity info
H_mean     = 18000000;
H_std      = 1000000;

Po_mean    = 1000;
Po_std     = 50;

F_s_mean   = 1400;
F_s_std    = 20;

theta_mean = 15;
theta_std  = 1;

L_sp_mean  = 2;
L_sp_std   = 0.4;

q_mean     = 0.5;
q_std      = 1;

RD_mean    = 5;
RD_std     = 1;

L_a_mean   = 2;
L_a_std    = 0.4;

C_d_mean   = 1;
C_d_std    = 0.3;

%% Setup Var_Info

var_info = struct;

% Save variable names in var_info
var_info.rand_inputs  = rand_inputs;
var_info.fixed_inputs = fixed_inputs;
var_info.output_vars  = output_vars;

% Save fixed variables
var_info.RE         = RE;
var_info.mu         = mu;
var_info.phi        = phi;
var_info.c          = c;
var_info.A_s        = A_s;
var_info.sun_i      = sun_i;
var_info.dt_slew    = dt_slew;
var_info.M          = M;
var_info.rho        = rho;
var_info.A          = A;
var_info.n          = n;
var_info.omega_max  = omega_max;
var_info.P_hold     = P_hold;
var_info.I_d        = I_d;
var_info.nu         = nu;
var_info.LT         = LT;
var_info.eps_deg    = eps_deg;
var_info.r_lw       = r_lw;
var_info.n_sa       = n_sa;
var_info.rho_sa     = rho_sa;
var_info.t          = t;
var_info.D          = D;
var_info.I_bodyx    = I_bodyx;
var_info.I_bodyy    = I_bodyy;
var_info.I_bodyz    = I_bodyz;

% Save mean and variance of random variables
var_info.H_mean      = H_mean;
var_info.Po_mean     = Po_mean;
var_info.F_s_mean    = F_s_mean;
var_info.theta_mean  = theta_mean;
var_info.L_sp_mean   = L_sp_mean;
var_info.q_mean      = q_mean;
var_info.RD_mean     = RD_mean;
var_info.L_a_mean    = L_a_mean;
var_info.C_d_mean    = C_d_mean;

var_info.H_std       = H_std;
var_info.Po_std      = Po_std;
var_info.F_s_std     = F_s_std;
var_info.theta_std   = theta_std;
var_info.L_sp_std    = L_sp_std;
var_info.q_std       = q_std;
var_info.RD_std      = RD_std;
var_info.L_a_std     = L_a_std;
var_info.C_d_std     = C_d_std;

% Save the mean outputs if they exist
if isfield(particle_system, 'ref_mean_out')

	% Setup coupling variables
	var_info.v 			= particle_system.ref_mean_out.v;
	var_info.dt_orbit 	= particle_system.ref_mean_out.dt_orbit;
	var_info.dt_eclipse = particle_system.ref_mean_out.dt_eclipse;
	var_info.PACS 		= particle_system.ref_mean_out.PACS;
	var_info.I_max		= particle_system.ref_mean_out.I_max;
	var_info.I_min 		= particle_system.ref_mean_out.I_min;
	var_info.dt_orbit   = particle_system.ref_mean_out.dt_orbit;
	var_info.theta_slew = particle_system.ref_mean_out.theta_slew;

end

end