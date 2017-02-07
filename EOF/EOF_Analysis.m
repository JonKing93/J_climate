function[s] = EOF_Analysis(Data, covcorr, MC, noiseType, confidence, varargin)
%% Performs a full EOF / PC Analysis of a data set.
% Finds the PC modes, explained variances, signals, and signals scaled to
% the standardized dataset. Performs a ruleN significance test on the modes
% and rotates the significant modes according to VARIMAX criterion.
%
% [s] = PCA_Analysis(Data, covcorr, MC, noiseType, confidence, analysisSpecs)
%
% ----- Inputs -----
%
% Data: A 2D data matrix. Each column corresponds to a particular time
%   series. Data cannot contain NaN entries.
%
% covcorr: The desired analysis matrix.
%       'cov': Covariance matrix
%       'corr': Correlation matrix
%       'none': Perform svd directly on data matrix. (The analysis will 
%               detrend but not zscore the data)
%  
% MC: The number of Monte Carlo iterations used in the ruleN significance
%       test
%
% noiseType: The noise used in the ruleN significance test
%       'white': white Gaussian noise
%       'red': lag-1 autocorrelated noise with added white noise
%
% confidence: The desired confidence interval. All modes passing the
%       confidence interval will be added to the VARIMAX rotation.
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
% s: A structure containing the following fields
%
%   Datax0: The standardized data matrix
%
%   C: The analysis matrix. This will be either the original data matrix or
%       its correlation or covariance matrix
%
%   eigVecs: The eigenvectors of the analysis matrix. These are the PCA
%       modes. Each column contains one mode.
%
%   eigVals: The eigenvalues of the analysis matrix. These are the loadings
%       of the PCA modes.
%
%   expVar: The explained variance of each PCA mode.
%
%   signals: The signals arising from each PCA mode. Each column is a
%       signal.
%
%   scaledSignals: The signals scaled to the standardized data matrix.
%
%   numSig: The number of PCA modes that pass the ruleN significance test.
%
%   randEigvals: The set of random, normalized eigenvalues generated during
%   the ruleN significance test. Each row contains the set of eigenvalues
%   at a particular confidence interval.
%
%   normEigvals: The normalized data eigenvalues
%
%   thresh: The index of the threshold row in randEigvals that data
%       eigenvalues must exceed to pass the ruleN significance test.
%
%   trueConf: The true confidence level of the threshold row.
%
%   scaledVecs: The significant eigenvectors scaled by the square root of
%       the eigenvalues. These are the values used for the VARIMAX
%       rotation.
%
%   rotatedModes: The VARIMAX rotated PCA modes.
%
%   rotatedEigvals: The eigenvalues / loadings for the rotated modes.
%
%   rotatedExpVar: The explained variance of the rotated eigenvalues.
%
%   rotationMatrix: The rotation matrix used to create the rotated modes.
%
%   rotatedSignals: The signals corresponding to the rotated modes
%
%   scaledRotSignals: The scaled signals from the rotated modes.


% Initial error checks
if nargin == 1
    extraCheck(Data)
    return;
end
[notFullSvd] = setup(covcorr, varargin{:});

% Declare the intial structure
s = struct();

% Run the PCA on the Data
[s.eigVals, s.eigVecs, s.Datax0, s.C] = simpleEOF(Data, covcorr, varargin{:});

% Get the signals
s.signals = getSignals(s.Datax0, s.eigVecs);

% Get the scaled signals
s.scaledSignals = scaleSignals(s.signals, s.eigVals);

% Get the explained variance
if notFullSvd   % Use the data variance if svd was incomplete
    s.expVar = explainedVar(s.eigVals, s.Datax0);
else
    s.expVar = explainedVar(s.eigVals);
end    

% Run the ruleN significance test
[s.numSig, s.randEigvals, s.normEigvals, s.thresh, s.conf] = ...
    ruleN(Data, s.eigVals, MC, noiseType, confidence, covcorr, varargin{:});

% Rotate the significant vectors (if more than 1 are significant)
if s.numSig < 2   % Less than 2 are significant, set rotation output to NaN
    s.scaledVecs = NaN;
    s.rotatedModes = NaN;
    s.rotationMatrix = NaN;
    s.rotatedEigvals = NaN;
    s.rotatedSignals = NaN;
    s.scaledRotSignals = NaN;
    
else % There are significant modes, rotate them...
    
    % Scale the eigenvectors
    s.scaledVecs = scaleEigvecs( s.eigVecs(:,1:s.numSig), s.eigVals(1:s.numSig));
    
    % Rotate the eigenvectors
    [s.rotatedModes, s.rotatedEigvals, s.rotatedExpVar, s.rotationMatrix] = ...
        varimaxRotation( s.scaledVecs, s.eigVals(1:s.numSig));
    
    % Get the rotated signals
    s.rotatedSignals = getSignals(s.Datax0, s.rotatedModes);
    
    % Get the scaled rotated signals
    s.scaledRotSignals = scaleSignals( s.rotatedSignals, s.rotatedEigvals);
    
end

end

%%%%% Helper Functions %%%%%
function[notFullSvd] = setup(covcorr, varargin)

notFullSvd = false;

% Determine if explainedVar will need extra inputs
for ks = 1:length(varargin)
    spec = varargin{ks};
    
    if isscalar(spec)  % Only the first several eigs are stored
        notFullSvd = true; 
    elseif strcmp(spec, 'econ') % Economy size SVD is performed
        notFullSvd = true;
    end
    
    if notFullSVD && strcmp(covcorr, 'none')
        error('Explained variance cannot be calculated when an incomplete SVD is performed directly on a data matrix.');
    end
end
    
% SVDS is always economy sized
if strcmp(covcorr, 'svds')
    notFullSvd = true;
end
    


end

function[] = extraCheck(egg)
% ... Happy Easter!
if strcmp(egg, 'Easter Egg') || strcmp(egg, 'easter egg') || strcmp(egg, 'Easter egg')
    fprintf('After careful analysis, PCs are better than Macs...\r\n');
else
    error('Insufficient input arguments for full PC Analysis');
end
end