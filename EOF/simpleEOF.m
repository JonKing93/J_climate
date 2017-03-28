function[modes, loadings, Datax0, C] = simpleEOF(Data, matrix, varargin)
%% Gets the EOF modes and loadings of a data matrix.
% 
% [loadings, modes, Datax0, C] = simpleEOF(Data, matrix)
%
% [...] = simpleEOF(Data, matrix, 'svds', 'econ')
% Uses an economy sized svds decomposition rather than the full svd.
%
% [...] = simpleEOF(Data, matrix, 'svds', nEigs)
% Uses the svds decomposition and determines the first nEigs eigenvalues.
%
% [...] = simpleEOF(Data, matrix, 'svd')
% Performs the normal svd.
%
%
% ----- Inputs -----
% 
% Data: A 2D data matrix. Each column corresponds to a particular time
%   series. Data cannot contain NaN entries.
%
% matrix: The desired analysis matrix.
%       'cov': Covariance matrix -- Minimizes variance along EOFs
%       'corr': Correlation matrix -- Minimizes relative variance along
%               EOFs. Often useful for data series with significantly
%               differenct magnitudes.
%       'none': Perform svd directly on data matrix. (This analysis will 
%               detrend but not zscore the data)
%
%
% ----- Outputs -----
%
% modes: The EOF modes. These are the eigenvectors of the analysis 
%       matrix. Each column is one mode.
%
% loadings: A vector of the mode loadings sorted from most to least 
%       significant. Loadings are the eigenvalues for the analysis matrix.
%
% modes: The EOF modes. These are the eigenvectors of the analysis 
%       matrix. Each column is one mode.
%
% Datax0: The standardized or detrended data matrix
%
% C: The analysis matrix for the PCA. This may be the stadardized data
%       matrix, its covariance matrix, or its correlation matrix, as appropriate.

% Get the inputs
[svdArgs] = parseInputs(varargin(:));

% Error check, determine analysis specifications
[covcorr, svdFunc, decompArg] = setup(Data, varargin{:});

% Standardize Data
Datax0 = standardizeData(Data, matrix);

% Get covariance OR correlation matrix OR leave data as is
C = getAnalysisMatrix( Datax0, matrix);

% Run SVD(S)
[loadings, modes] = quickSVD(C, svdFunc, decompArg);

% Calculate eigenvalues if SVD is performed directly on data
if strcmp(covcorr, 'none')
    loadings = (loadings.^2) / (length(loadings) - 1);
end

end


%%%%% Helper Functions %%%%%
function[svdArgs] = parseInputs(inArgs)

% Set the default
svdArgs = {'svd'};

if ~isempty( inArgs)
    for k = 1:length(inArgs)
        arg = inArgs{k};
        
        if strcmpi(arg, 'svd')
            % Do nothing
        elseif strcmpi(arg, 'svds')
            if length(inArgs) >= k+1 && ( isscalar(inArgs{k+1}) || strcmpi(inArgs{k+1},'econ') )
                svdArgs = {'svds', inArgs{k+1}};
            else
                error('The svds flag must be followed by nEigs or the ''econ'' flag');
            end
        else
            error('Unrecognized Input');
        end
    end
end
end


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