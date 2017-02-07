function[eigVals, eigVecs, Datax0, C] = simpleEOF(Data, varargin)
%% Gets the PCs (EOFs) of a data matrix
%
% [eigVals, eigVecs, Datax0, C] = simplePCA(Data, analysisSpecs)      
%
% ----- Inputs -----
% 
% Data: A 2D data matrix. Each column corresponds to a particular time
%   series. Data cannot contain NaN entries.
%
% *** Optional Inputs ***
%
% analysisSpecs: Additional specifications for the covariance / correlation
%       matrices, and the svd decomposition
%
%   Covariance / Correlation Matrices:
%           'cov': (default) Uses the covariance matrix
%           'corr': Uses the correlation matrix
%           'none': Will perform svd directly on data matrix. (The analysis will 
%                   detrend but not zscore the data)
%        
%   Decomposition type: svd or svds (these may have differing runtimes)
%           'svd': (default) Uses the svd function for decomposition
%           'svds': Uses the svds function for decomposition
%             
%   Decomposition size: The number of eigenvectors found (may affect runtime)
%           'all': (default)Performs the full decomposition
%           'econ': Performs the economy size decomposition
%           neigs: An integer specifying the number of leading eigenvectors to find                
%
%
% ----- Outputs -----
%
% eigVals: a vector of the desired eigenvalues sorted from most to least 
%     significant.
%
% eigVecs: The eigenvectors of the corr/cov matrix
%
% Datax0: Standardized data
%    mean = 0 if the covariance matrix is used
%    variance = 1, mean = 0 if the correlation matrix is used
%
% C: The analysis matrix for the PCA. This may be the stadardized data
%   matrix, its covariance matrix, or its correlation matrix, as appropriate.
%
% ----- Suggested Readings -----
% 
% IN PROGRESS...
%
%

% Error check, determine analysis specifications
[covcorr, svdFunc, decompArg] = setup(Data, varargin{:});

% Standardize Data
Datax0 = standardizeData(Data,covcorr);

% Get covariance OR correlation matrix OR leave data as is
C = getAnalysisMatrix( Datax0, covcorr);

% Run SVD(S)
[eigVals, eigVecs] = quickSVD(C, svdFunc, decompArg);

% Calculate eigenvalues if SVD is performed directly on data
if strcmp(covcorr, 'none')
    eigVals = (eigVals.^2) / (length(eigVals) - 1);
end

end


%%%%% Helper Functions %%%%%
function[covcorr, svdFunc, decompArg] = setup(Data, varargin)
%% Setup function for simple PCA

%% Ensure data matrix is 2D    
if ~ismatrix(Data)
    error('Data must be a 2D matrix');
end

% Ensure data does not contain NaNs
if NaNcheck(Data)
    error('Data cannot contain NaNs');
end

% Set defaults
covcorr = 'cov';
svdFunc = 'svd';
decompArg = 'all';

% Return if no additional specifications
if isempty(varargin)
    return;
end

%% Determine what to do with different analysis specifications
for k = 1:length(varargin)
    spec = varargin{:};  
    
    if strcmp( spec, 'cov') % Use covariance for cov string
        covcorr = 'cov';
    
    elseif strcmp( spec, 'corr') % Use correlation for corr string
        covcorr = 'corr';
        
    elseif strcmp( spec, 'none') % Use no cov / corr matrix
        covcorr = 'none';
        
    elseif strcmp( spec, 'svd' ) % Use svd
        svdFunc = 'svd';
        
    elseif strcmp( spec, 'svds' ) % Use svds
        svdFunc = 'svds';
        
    elseif strcmp( spec, 'all' ) % Full svd
        decompArg = 'all';
        
    elseif strcmp( spec, 'econ' ) % Economy svd
        decompArg = 'econ';
        
    elseif isscalar(spec)    % Specific number of eigs
        decompArg = spec;
        
    else % Error for anything else
        error('Unrecognized Input');
    end    
end

end

function[Datax0] = standardizeData(Data, covcorr)
% Standardizes the data as appropriate
if strcmp(covcorr, 'corr')
    Datax0 = zscore(Data);
else
    Datax0 = detrend(Data, 'constant');
end
end

function[C] = getAnalysisMatrix( Datax0, covcorr)
if strcmp(covcorr, 'cov')
    C = cov(Datax0);
    
elseif  strcmp(covcorr,'corr')
    C = corr(Datax0);
else
    C = Datax0;
end
end