function [F] = satellite_solver(outputs, var_info, coupling_vars)
% SATELLITE_SOLVER: Solves satellite problem using fixed point iteration
% for given det_inputs and rand_inputs.
%
% See COPYRIGHT for license information 
%
% The satellite problem is adapted from the reference:
% Sankararaman, S., and Mahadevan, S. Likelihood-based approach to multidisciplinary
% analysis under uncertainty. Journal of Mechanical Design, 2012.

% Setup fixed input variables
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

% Setup random inputs variables
H          = var_info.H;
Po         = var_info.Po;
F_s        = var_info.F_s;
theta      = var_info.theta;
L_sp       = var_info.L_sp;
q          = var_info.q;
RD         = var_info.RD;
L_a        = var_info.L_a;
C_d        = var_info.C_d;

% Setup orbit output variables
v          = outputs(1);
dt_orbit   = outputs(2);
dt_eclipse = outputs(3);
theta_slew = outputs(4);

% Set attitude control output variables
tau_slew   = outputs(5);
tau_g      = outputs(6);
tau_sp     = outputs(7);
tau_m      = outputs(8);
tau_a      = outputs(9);
tau_dist   = outputs(10);
tau_tot    = outputs(11);
PACS       = outputs(12);

% Set power output variables
P_tot      = outputs(13);
T_e        = outputs(14);
T_d        = outputs(15);
P_e        = outputs(16);
P_d        = outputs(17);
P_sa       = outputs(18);
PBOL       = outputs(19);
PEOL       = outputs(20);
A_sa       = outputs(21);
L          = outputs(22);
W          = outputs(23);
m_sa       = outputs(24);
I_sax      = outputs(25);
I_say      = outputs(26);
I_saz      = outputs(27);
I_max      = outputs(28);
I_min      = outputs(29);

%% Compute Orbit Subsystem Outputs

% Setup inputs
orbit = struct;
orbit.mu  = mu;
orbit.RE  = RE;
orbit.H   = H;
orbit.phi = phi;

% Compute output quantitites
orbit.v          = sqrt(orbit.mu/(orbit.RE + orbit.H));
orbit.dt_orbit   = 2*pi*(orbit.RE + orbit.H)/v;
orbit.dt_eclipse = orbit.dt_orbit/pi*asin(orbit.RE/(orbit.RE+orbit.H));
orbit.theta_slew = atan(sin(orbit.phi/orbit.RE)/(1 - cos(orbit.phi/orbit.RE) + orbit.H/orbit.RE));

% Assemble Orbit outputs
orbit_out  = [orbit.v, orbit.dt_orbit, orbit.dt_eclipse, orbit.theta_slew];

%% Compute Attitude Control Subsystem Outputs

% Setup inputs
ac = struct;
ac.dt_slew   = dt_slew;
ac.mu        = mu;
ac.sun_i     = sun_i;
ac.P_hold    = P_hold;
ac.n         = n;
ac.omega_max = omega_max;
ac.L_a       = L_a;
ac.rho       = rho;
ac.C_d       = C_d;
ac.A         = A;
ac.M         = M;
ac.RD        = RD;
ac.RE        = RE;
ac.H         = H;
ac.L_sp      = L_sp;
ac.F_s       = F_s;
ac.c         = c;
ac.A_s       = A_s;
ac.q         = q;
ac.theta     = theta;

% Setup coupling variables
if coupling_vars(4) == 1
    ac.I_max = I_max;
else
    ac.I_max = var_info.I_max;
end

if coupling_vars(5) == 1
    ac.I_min = I_min;
else
    ac.I_min = var_info.I_min;
end

if coupling_vars(6) == 1
    ac.dt_orbit = orbit.dt_orbit;
else
    ac.dt_orbit = var_info.dt_orbit;
end

if coupling_vars(7) == 1
    ac.v = orbit.v;
else
    ac.v = var_info.v;
end

if coupling_vars(8) == 1
    ac.theta_slew = orbit.theta_slew;
else
    ac.theta_slew = var_info.theta_slew;
end

% Compute output quantitites
ac.tau_slew   = 4*ac.theta_slew/(ac.dt_slew)^2*ac.I_max;
ac.tau_g      = 3*ac.mu/(2*(ac.RE+ac.H)^3)*abs(ac.I_max - ac.I_min)*sin(2*ac.theta);
ac.tau_sp     = ac.L_sp*ac.F_s/ac.c*ac.A_s*(1+ac.q)*cos(ac.sun_i);
ac.tau_m      = 2*ac.M*ac.RD/(ac.RE+ac.H)^3;
ac.tau_a      = 0.5*ac.L_a*ac.rho*ac.C_d*ac.A*ac.v^2;
ac.tau_dist   = sqrt(ac.tau_g^2 + ac.tau_sp^2 + ac.tau_m^2 + ac.tau_a^2);
ac.tau_tot    = max(ac.tau_slew,ac.tau_dist);
ac.PACS       = ac.tau_tot*ac.omega_max + ac.n*ac.P_hold;

% Assemble AC outputs
ac_out = [ac.tau_slew, ac.tau_g, ac.tau_sp, ac.tau_m, ac.tau_a, ...
          ac.tau_dist, ac.tau_tot, ac.PACS];

%% Compute Power Subsystem Outputs

% Setup external inputs
power = struct;
power.I_bodyx = I_bodyx;
power.I_bodyy = I_bodyy;
power.I_bodyz = I_bodyz;
power.D       = D;
power.t       = t;
power.rho_sa  = rho_sa;
power.A_sa    = A_sa;
power.r_lw    = r_lw;
power.n_sa    = n_sa;
power.eps_deg = eps_deg;
power.LT      = LT;
power.F_s     = F_s;
power.I_d     = I_d;
power.sun_i   = sun_i;
power.Po      = Po;
power.nu      = nu;

% Setup coupling variables
if coupling_vars(1) == 1
    power.dt_orbit = orbit.dt_orbit;
else
    power.dt_orbit = var_info.dt_orbit;
end

if coupling_vars(2) == 1
    power.dt_eclipse = orbit.dt_eclipse;
else
    power.dt_eclipse = var_info.dt_eclipse;
end

if coupling_vars(3) == 1
    power.PACS = ac.PACS;
else
    power.PACS = var_info.PACS;
end

% Compute output quantitites
power.P_tot  = power.PACS + power.Po;
power.T_e    = power.dt_eclipse;
power.T_d    = power.dt_orbit - power.T_e;
power.P_e    = power.P_tot;
power.P_d    = power.P_tot;
power.P_sa   = (power.P_e*power.T_e/0.6 + power.P_d*power.T_d/0.8)/power.T_d;
power.PBOL   = power.nu*power.F_s*power.I_d*cos(power.sun_i);
power.PEOL   = power.PBOL*(1-power.eps_deg)^(power.LT);
power.A_sa   = power.P_sa/power.PEOL;

% Compute solar array geometry quantities
power.L      = sqrt(power.A_sa*power.r_lw/power.n_sa);
power.W      = sqrt(power.A_sa/(power.r_lw*power.n_sa));
power.m_sa   = 2*power.rho_sa*power.L*power.W*power.t;

% Compute solar array moments of inertia
power.I_sax  = power.m_sa*(1/12*(power.L^2 + power.t^2) + (power.D + power.L/2)^2);
power.I_say  = power.m_sa/12*(power.t^2 + power.W^2);
power.I_saz  = power.m_sa*(1/12*(power.L^2 + power.W^2) + (power.D + power.L/2)^2);
power.I_tot  = [power.I_sax, power.I_say, power.I_saz] + [power.I_bodyx, power.I_bodyy, power.I_bodyz];
power.I_max  = max(power.I_tot);
power.I_min  = min(power.I_tot);
 
% Assemble power outputs
power_out    = [power.P_tot, power.T_e, power.T_d, power.P_e, power.P_d, ...
                power.P_sa, power.PBOL, power.PEOL, power.A_sa, power.L, ...
                power.W, power.m_sa, power.I_sax, power.I_say, power.I_saz, ...
                power.I_max, power.I_min];

%% Set Outputs

% Set output quantitites
F = outputs - [orbit_out, ac_out, power_out]';

% -- END OF FILE -- 