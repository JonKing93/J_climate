% EOF
%
% This folder contains functions involved in EOF Analysis.
%
% The function EOF_Analysis runs an analysis in full, but the functions may
% also be used individually.
%
% Files
%   PC_Analysis  - Performs a full PCA (EOF) Analysis of a data set.
%   simplePCA    - Gets the PCs (EOFs) of a data matrix
%   explainedVar - Gets the percent explained variance of a set of eigenvalues / PC loadings
%   getSignals   - Gets the signals from a standardized dataset and its eigenvectors.
%   scaleSignals - Scales signals to the standard deviation of the standardized data.
%   ruleN        - Runs a Rule N significance test on a data matrix and its eigenvalues / PC loadings.
%   scaleEigvecs - Scales the eigenvectors by the square root of the eigenvalues.
%   varimaxRotation - Performs a VARIMAX rotation on a set of scaled eigenvectors and eigenvalues.
%
% Written by: Jonathan King (jonking93@email.arizona.edu)
% V1.0  (11-17-2016)
%
% This work based on the course Spatiotemporal Data Analysis, presented by
% Kevin Anchukaitis.
