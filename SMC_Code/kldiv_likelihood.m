function kl = kldiv_likelihood(particle_system, model)
% KLDIV_LIKELIHOOD: Function computes the KL divergence between the
% linearized version of the models to use as an approximation for the
% likelihood term in the SMC algorithm
%
% See COPYRIGHT for license information 
%

% Extract inputs from particle_system
ref_mean   = particle_system.ref_mean;
ref_cov    = particle_system.ref_cov;
n_vars     = particle_system.n_vars;

%% Run Gaussian Linearization

% Determine the dimension of QoI
n_dim = length(ref_mean);

% Open model_data.txt file
fileID = fopen('Model/model_data.txt','r');
model_data = textscan(fileID, repmat('%f ', 1, n_vars+n_dim+n_dim^2));
model_data = cell2mat(model_data);
fclose(fileID);

% Extract existing models
model_runs = model_data(:,1:n_vars);
[model_exists, model_idx] = ismember(model, model_runs, 'rows');

% Check if model has been evaluated in model_dict
if (model_exists==1)

    % Extract mean and covariance
    dec_mean = model_data(model_idx,n_vars+1:n_vars+n_dim);
    dec_cov  = model_data(model_idx,n_vars+n_dim+1:end);

    % Reshape mean vector and covariance matrix
    dec_mean = dec_mean';
    dec_cov  = reshape(dec_cov, n_dim, n_dim);

else

    % Find mean and covariance of decoupled model
    [dec_mean, dec_cov] = gaussian_linearization(particle_system, model);
    
    % Store result in model_data file
    fileID = fopen('Model/model_data.txt','a');
    fprintf(fileID, '%.16e ', [model, dec_mean', reshape(dec_cov, 1, n_dim^2)]);
    fprintf(fileID, '\n');
    fclose(fileID);

end

%% Compute KL Divergence

% Compute KL divergence between Gaussians
kl = kldiv_mvn(ref_mean, dec_mean, ref_cov, dec_cov);

% If Kl is < 0, set to Infinity
if (kl < 0) || (imag(kl) ~= 0) || isnan(kl)
    kl = Inf;
end

end
