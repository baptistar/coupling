% RUN_MODEL_SELECT: Script runs the model selection algorithm
%
% See COPYRIGHT for license information 
%

% Generate file to store model_data
fclose(fopen('Model/model_data.txt', 'w'));

% Generate folder to store Results
mkdir('Results')

% Find the number of lambda values
n_lambda = length(lambda_vect);

% Declare cell to store optimal model
opt_models = cell(n_lambda,1);

% Run model selection algorithm for each lambda
for i=1:n_lambda

	% Setup input_struct
	input_struct = struct;
	input_struct.eng_problem = eng_problem;
	input_struct.n_models  = n_particles;
	input_struct.delta     = particle_div;
	input_struct.lambda    = lambda_vect(i);
	input_struct.n_vars    = n_coupling;
	input_struct.qoi       = prob_qoi;
	
	% Find optimal model
	opt_models{i} = binary_SMC(input_struct);

end
