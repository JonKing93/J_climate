function[lastSigNum, randEigSort, normEigvals, thresh, realConf] = ...
    ruleN(Data, eigVals, MC, noiseType, confidence, covcorr, varargin)
%% Runs a Rule N significance test on a data matrix and its eigenvalues / PC loadings.
%
% [lastSigNum, randEigSort, normEigvals, thresh, realConf] = ...
%    ruleN(Data, eigVals, MC, noiseType, confidence, covcorr, analysisSpecs)
%
% ----- Inputs -----
%
% Data: A 2D data matrix. Each column corresponds to a series of
%   observations.
%
% normVals: The eigenvalues of Data
%
% MC: The Monte Carlo number of iterations
%
% noiseType: 
%   'white':    white noise
%   'red':      lag-1 red noise with added white noise
%
% confidence: Desired confidence interval for ruleN to pass. Must be on the
%       interval (0 1)
%
% covcorr: The desired analysis matrix.
%       'cov': Covariance matrix
%       'corr': Correlation matrix
%       'none': Perform svd directly on data matrix. (The analysis will 
%               detrend but not zscore the data)
%
% *** Optional Inputs ***
%
% analysisSpecs:
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
% lastSigNum: The number of eigenvalues that pass rule N
%
% randEigSort: The matrix of random, normalized, sorted eigenvalues
%
% normEigvals: The normalized data eigenvalues
% 
% thresh: The integer threshold that eigenvalues were required to pass
%
% realConf: The true confidence interval of this threshold
% 
%
% ----- Background Reading -----
%
% Principal Component Analysis in Metereology and Oceanography. Rudolph
% Preisendorfer. Elsevier Science Publishers. New York, 1988.
%

[ar1, normEigvals] = setup(Data, eigVals, noiseType, confidence, MC);

% Get the data size
[m, n] = size(Data);

% Run rule N
randEigvals = NaN(MC,n);
for k = 1:MC
    
    % Create a random matrix
    g = buildMatrix(noiseType,m,n,ar1);
    
    % Scale to the standard deviation of the original matrix
    g = g * sqrt( diag( var( Data)));
    
    % Run an EOF analysis on the random matrix
    [randEig, ~] = simplePCA(g, covcorr);
    
    % Normalize the eigenvalues
    randEig = randEig ./ sum(randEig);
    
    % Store the random eigenvalues
    randEigvals(k,:) = randEig;
    
end

% Sort the eigenvalues
randEigSort = sort(randEigvals);

% Calculate the confidence level threshold
thresh = ceil( MC * confidence);
realConf = MC / thresh;

% Find the significant values
for k = 1:n
    if normEigvals(k) <= randEigSort(thresh, k)
        lastSigNum = k-1;
        break;
    end
end

end

%%%%% Helper Functions %%%%%
function[g] = buildMatrix(noiseType, m, n, ar1)
%% Builds the matrix g as appropriate for red or white noise
switch noiseType
    
    % Random matrix for white noise
    case 'white'
        % Create a random matrix
        g = randn(m,n);
    
    % Add lag-1 autocorrelation for red noise
    case 'red'
        % Preallocate 
        g = NaN(m,n);
        
        % Create random first row
        g(1,:) = randn(1,n);
        
        % Calculate autocorrelation through matrix. Add random noise
        for j = 1:m-1
            g(j+1,:) = (ar1' .* g(j,:)) + randn(1,n);
        end
        
        % Standardize so later scaling is correct
        g = zscore(g);   
end
end

function[ar1, normEigvals] = setup(Data, eigVals, noiseType, confidence, MC)

% Ensure Data is 2D
if ~ismatrix(Data)
    error('RuleN is for 2D Data matrices');
end

% Get noise type
if ~( strcmp(noiseType,'red') || strcmp(noiseType, 'white') )
    error('Unrecognized noise type');
end

% Precalculate ar1 if required
if strcmp(noiseType, 'red')
    r = corr( Data(1:end-1,:), Data(2:end,:) );
    ar1 = diag(r);
else
    ar1 = NaN;
end    

% Normalize Eigenvalues
normEigvals = eigVals ./ sum(eigVals);

% Ensure confidence interval is on (0 1)
if confidence <=0 || confidence >=1
    error('confidence must be on the interval (0,1)');
end

% Ensure the Monte Carlo number is positive
if MC < 1
    error('The Monte Carlo number must be a positive integer');
end


end
    