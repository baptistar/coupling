% MODEL_SELECTION
% This code implements automatic coupling model selection for
% multidisciplinary engineering problems
%
% For a reference to the algorithmic formulation see:
% Baptista, R., Willcox, K., Marzouk, Y., and Peherstorfer, B.,
% Optimal Approximations of Coupling in Multidisciplinary Models,
% AIAA Scitech Proceedings, January 2017
%
% Author:  Ricardo Santos-Baptista 
% Date:    November 2016
%
% See COPYRIGHT for license information 
%

clear; close all; clc

addpath SMC_Code
addpath Aux_Functions
addpath Model

% Setup algorithm parameters
eng_problem  = @satellite_model;
prob_qoi     = {'tau_tot','P_tot','A_sa'};
lambda_vect  = [1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4];
n_particles  = 100;
particle_div = 0.1;
n_coupling   = 8;

% Run model selection algorithm
run_model_select;

% -- END OF FILE --
