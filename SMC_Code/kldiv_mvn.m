function D_KL = kldiv_mvn(mean_1, mean_2, cov_1, cov_2)
% KLDIV_MVN: Fucntion computes the KL divergence between two multivariate
% Gaussian distributions as specified by their mean vectors and covariance
% matrices.
%
% See COPYRIGHT for license information 
%

warning('off','MATLAB:singularMatrix')
warning('off','MATLAB:nearlySingularMatrix')

% Determine the size of the vectors
K = length(mean_1);

% Compute KL divergence
D_KL = 0.5*(trace(cov_2\cov_1) + (mean_2 - mean_1)'*(cov_2\(mean_2 - mean_1)) ...
       - K + log(det(cov_2)/det(cov_1)));

end