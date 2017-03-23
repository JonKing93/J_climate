function[s] = EOF_Analysis(Data, matrix, MC, noiseType, pval, varargin)
%% Performs a full EOF Analysis of a data set.
% 
% [s] = EOF_Analysis(Data, matrix, MC, noiseType, pval)
% Finds the EOF modes, explained variances, signals, and signals scaled to
% the standardized dataset. Performs a Rule N significance test on the modes
% and rotates the significant modes according to VARIMAX criterion. Returns
% all calculated values in a structure, s.
%
% [s] = EOF_Analysis(..., 'testMC')
% Also saves the eigenvalues at the required significance level for each
% Monte Carlo iteration to test convergence.
%
% [s] = EOF_Analysis(..., 'showProgress')
% Displays the current Monte Carlo iteration onscreen.
%
% [s] = EOF_Analysis(..., 'svds', 'econ')
% Performs the economy sized svds decomposition rather than the default
% svd.
%
% [s] = EOF_Analysis(..., 'svds', nEigs)
% Uses the svds decomposition and determines the first nEigs eigenvalues.
%
%
% ----- Inputs -----
%
% Data: A 2D data matrix. Each column corresponds to a particular data
%       series. Data may only contain numeric entries.
%
% matrix: The desired analysis matrix.
%       'cov': Covariance matrix -- Minimizes variance along EOFs
%       'corr': Correlation matrix -- Minimizes relative variance along
%               EOFs. Often useful for data series with significantly
%               differenct magnitudes.
%       'none': Perform svd directly on data matrix. (This analysis will 
%               detrend but not zscore the data)
%  
% MC: The number of Monte Carlo iterations used in the Rule N significance
%       test
%
% noiseType: The noise used in the ruleN significance test
%       'white': white Gaussian noise
%       'red': lag-1 autocorrelated noise with added white noise
%
% pval: The significance level that the significance test should pass.        
%
%
% ----- Outputs -----
%
% s: A structure containing the following fields
%
%   Datax0: The standardized or detrended data matrix
%
%   C: The analysis matrix. The covariance or correlation matrix of Datax0,
%       or the original dataset.
%
%   modes: The EOF modes. These are the eigenvectors of the analysis 
%       matrix. Each column is one mode.
%
%   loadings: The loadings of each data series on each EOF modes. These are
%       the eigenvalues of the analysis matrix. Each column corresponds to
%       a different mode.
%
%   normLoads: The normalized loadings.
%
%   expVar: The data variance explained by each mode.
%
%   signals: The signal from each EOF mode. Signals are the imprint of each
%       mode on the original data series. Each column is a signal.
%
%   scaSignals: The signals scaled to the standardized data matrix. Allows
%       for quick comparison of signals with the original data series.
%
%   nSig: The number of modes that pass the Rule N significance test. 
%
%   randLoads: The set of random, normalized loadings generated during
%       the Rule N significance test. Each row contains the set of loadings
%       at a particular confidence interval.
%
%   thresh: The index of the threshold row of random loadings in randLoads
%       that the original data series loadings must exceed in order to pass
%       the significance test.
%
%   trueSig: The true significance level of the threshold row.
%
%   iterTrueSig: The true significance level of the threshold row after
%       each iteration of the Monte Carlo simulations.
%
%   iterSigLevel: The set of loadings for each iteration of the Monte Carlo
%       simulations that the data series loadings must exceed in order to
%       pass the significance test.
%
%   scaModes: The scaled modes used for VARIMAX rotation. Modes are scaled 
%       by the square root of the loadings.
%
%   rotModes: The VARIMAX rotated modes.
%
%   rotLoads: The loadings for the rotated modes.
%
%   rotExpVar: The variance explained by the rotated loadings.
%
%   rotSignals: The signals corresponding to the rotated modes.
%
%   scaRotSignals: The scaled signal for each rotated mode.
%
%   rotMatrix: The rotation matrix used to rotate the selected.
%
%   metadata: Information concerning the settings used for the analysis.

% Initial error checks


[notFullSvd] = setup(matrix, varargin{:});

% Declare the intial structure
s = struct();

% Run the initial EOF on the Data
[s.eigVals, s.eigVecs, s.Datax0, s.C] = simpleEOF(Data, matrix, varargin{:});

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
[s.numSig, s.randEigvals, s.normEigvals, s.thresh, s.conf, s.iterTrueConf, s.iterConfEigs] = ...
    ruleN(Data, s.eigVals, MC, noiseType, pval, matrix, varargin{:});

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

