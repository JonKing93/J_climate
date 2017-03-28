function[eigVals, modes, Datax0, C] = simpleEOF(Data, matrix, varargin)
%% Gets the EOF modes and loadings of a data matrix.
% 
% [eigVals, modes, Datax0, C] = simpleEOF(Data, matrix)
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
%               EOFs. Useful for data series with significantly
%               different magnitudes.
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
% C: The analysis matrix for the PCA. This may be the standardized data
%       matrix, its covariance matrix, or its correlation matrix, as appropriate.

% Get the inputs
[svdArgs] = parseInputs(varargin(:));

% Error check
errCheck(Data, matrix);

% Standardize Data
Datax0 = standardizeData(Data, matrix); % Helper function

% Get the covariance / correlation / detrended data matrix.
C = getAnalysisMatrix( Datax0, matrix);

% Run SVD(S)
[eigVals, modes] = quickSVD(C, svdArgs);

% Calculate eigenvalues if SVD is performed directly on data
if strcmpi(matrix, 'none')
    eigVals = (eigVals.^2) / (length(eigVals) - 1);
end

end


%%%%% Helper Functions %%%%%
function[svdArgs] = parseInputs(varargin)
inArgs = varargin;

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

function[] = errCheck(Data, matrix)
%% Ensure data matrix is 2D    
if ~ismatrix(Data)
    error('Data must be a 2D matrix');
end

% Ensure data does not contain NaNs
if NaNcheck(Data)
    error('Data cannot contain NaNs');
end

% Matrix is recognized
if ~any( strcmpi(matrix, {'corr','cov','none'}) )
    error('Unrecognized matrix');
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

function[C] = getAnalysisMatrix( Datax0, matrix)
if strcmpi(matrix, 'cov') || strcmpi(matrix, 'corr')
    % Note that the corr matrix has been zscored, so taking the covariance
    % is equivalent to the correlation matrix.
    C = cov(Datax0);
else
    C = Datax0;
end
end