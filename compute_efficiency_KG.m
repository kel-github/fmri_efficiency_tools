function eff = compute_efficiency_KG( design_mat, contrast, sessions, nHRFB )

% design mat = full design matrix (e.g. from spm )
% contrast   = contrast matrix
% sessions   = number of sessions/runs
% nHRFB      = number of basis functions for HRF

contrast  = repmat( contrast, 1, sessions ); % across sessions
contrast   = kron(contrast,eye(nHRFB)/nHRFB); % for each HRF parameter
contrast   = cat(2, contrast, zeros(size(contrast, 1), 2));
iXX = pinv(design_mat'*design_mat); % compute pseudo inverse (covariance)
eff = 1/trace(contrast*iXX*contrast'); % compute efficiency





