function[expVar] = explainedVar(eigVals, varargin)
%% Gets the percent explained variance of a set of eigenvalues / PC loadings
%
% [expVar] = explainedVar(eigVals, Datax0)
%
% ----- Inputs -----
%
% eigVals: A set eigenvalues. If eigVals is a matrix, each column contains 
%       a set of eigenvalues.
%
%       Note: eigVals is assumed to contain all of the singular values for
%       a matrix. See the optional inputs section if this is not the case.
%
% *** Optional Inputs ***
%
% Datax0: The standardized matrix for the set of eigenvalues. If Datax0 is
%       3D, each dim1 x dim2 matrix is a standardized data matrix.
%
%       For economy svd or selection of just a few eigs in an svds, the 
%       reduced set of eigenvectors is insufficient to calculate explained
%       variance. However, mathematically, the eigenvalues must contain the
%       full amount of original data variance, albeit redistributed. Thus,
%       knowledge of total data variance can allow calculation of explained
%       variance. 

%       However, this DOES NOT APPLY to direct svd of a data matrix. A full
%       svd should be performed when directly analyzing data matrices to
%       allow calculation of explained variance.
%
%       Note: If a correlation matrix is used in PCA, the total variance of
%       the original data matrix will not give correct explained variance.
%       Thus, the standardized data matrix should be used.
%
%
% ----- Outputs -----
%
% expVar: a matrix containing the percent explained variance of each 
%   eigenvalue. Each column contains the explained variances of one set of
%   eigenvalues in the eigVals matrix.
%

% Run an error check, get the total variance from standardized data
[dataVar, neigs] = setup(eigVals, varargin{:});

% Get the appropriate total variance
if isnan(dataVar) % Sum of eigs for complete set of eigenvalues
    totalVar = sum(eigVals,1);
else                % Sum of data variance for 
    totalVar = dataVar;
end

expVar = 100 * eigVals ./ repmat(totalVar, neigs,1);

end

%%%%% Helper Functions %%%%%
function[dataVar, nEigs] = setup(eigVals, varargin)

% Ensure the right number of input args is given
if nargin > 2
    error('Too many input arguments');
end

% Ensure eigVals is a matrix
if ~ismatrix(eigVals)
    error('eigVals must be a matrix');
end

% Check for NaNs
if NaNcheck(eigVals)
    error('eigVals should not contain NaN');
end

% If dataVar is input...
if ~isempty(varargin)
    Datax0 = varargin{1};
    
    % Ensure dataVar is 3D
    if ndims(Datax0) > 3
        error('Datax0 must be 3D');
    end
    
    nSets = size(Datax0,3);    
    % Ensure dataVar has the correct size
    if nSets ~= size(eigVals,2)
        error('The number of entries in Datax0 does not match the number of sets of eigenvalues');
    end
    
    % Check for NaN
    if NaNcheck(Datax0)
        error('Datax0 cannot contain NaN');
    end
    
    % Calculate data variance
    dataVar = NaN(nSets,1);
    for k = 1:nSets
        dataVar(k) = sum(var(Datax0(:,:,k)));
    end        
    
    % Check that the sum of sets of eigenvalues do not exceed the values in
    % dataVar
    if any( sum(eigVals, 1) > dataVar )
        error('A set of eigenvalues exceeds the total variance within the variance of Datax0.\r\n');
    end
else
    dataVar = NaN;
end

% Get the number of eigenvectors per column
nEigs = size(eigVals,1);
end 
    

