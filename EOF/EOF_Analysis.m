function[s] = EOF_Analysis(Data, matrix, MC, noiseType, pval, varargin)
%% Performs a full EOF Analysis of a data set.
% 
% [s] = EOF_Analysis(Data, matrix, MC, noiseType, pval)
% Finds the EOF modes, explained variances, signals, and signals scaled to
% the standardized dataset. Performs a Rule N significance test on the modes
% and rotates the significant modes according to VARIMAX criterion. Returns
% all calculated values in a structure, s.
%
% [s] = EOF_Analysis(..., 'showProgress')
% Displays the current Monte Carlo iteration onscreen.
%
% [s] = EOF_Analysis(..., 'svds', 'econ')
% Performs the economy sized svds decomposition rather than the default
% svd.
%
% [s] = EOF_Analysis(..., 'svds', nModes)
% Uses the svds decomposition and determines the first nModes modes.
%
% [s] = EOF_Analysis(..., 'noSigTest')
% A flag to block the Rule N significance testing. The returned structure
% will not contain any fields requiring the significance test.
%
% [s] = EOF_Analysis(..., 'noConvergeTest')
% A flag to block the test for Monte Carlo convergence. The returned
% structure will contain neither the iterSigLevel nor iterTrueSig fields.
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
%               different magnitudes.
%       'none': Perform svd directly on data matrix.
%  
% MC: The number of Monte Carlo iterations used in the Rule N significance
%       test
%
% noiseType: The noise used in the Rule N significance test
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
%   eigVals: A vector with the eigenvalues of the analysis matrix. Each 
%       eigenvalue corresponds to an EOF mode. The larger the eigenvalue
%       magnitude, the more data variance explained by the mode.
%
%   modes: The EOF modes. These are the eigenvectors of the analysis 
%       matrix, also known as loadings. Each column is one mode.
%
%   expVar: The data variance explained by each mode.
%
%   signals: The signal for each EOF mode. Signals are the imprint of each
%       mode on the original data series, also known as scores or EOF 
%       time series. Each column is one signal.
%
%   scaSignals: The signals scaled to the standardized data matrix. Allows
%       for quick comparison of signals with the original data series.
%
%   normVals: The normalized data eigenvalues.
%
%   nSig: The number of modes that pass the Rule N significance test. 
%
%   randEigvals: The set of random, normalized eigenvalues generated during
%       the Rule N significance test. Each row contains the eigenvalues
%       at a particular confidence interval.
%
%   thresh: The index of the threshold row in randEigvals that the data 
%       eigenvalues must exceed in order to pass the significance test.
%
%   trueConf: The true confidence level of the threshold row.
%
%   iterTrueConf: The true confidence level of the threshold row after
%       each iteration of the Monte Carlo simulations.
%
%   iterSigEigs: The set of eigenvalues that the data values must exceed
%       for significance after each successive Monte Carlo iteration.
%
%   scaModes: The scaled modes used for VARIMAX rotation. Modes are scaled 
%       by the square root of the loadings.
%
%   rotModes: The VARIMAX rotated modes.
%
%   rotEigvals: The eigenvalues for the rotated modes.
%
%   rotExpVar: The variance explained by the rotated loadings.
%
%   rotSignals: The signals corresponding to the rotated modes.
%
%   scaRotSignals: The scaled signal for each rotated mode.
%
%   rotMatrix: The rotation matrix used to rotate the significant modes.
%
%   metadata: Information concerning the settings used for the analysis.
%       Contains: matrix, MC, noisetype, pval, and any additional flags.

% Parse Inputs, all error checking will occur in called functions
[showProgress, svdArgs, blockMC, convergeTest] = parseInputs(varargin(:));

% Do some setup
[incompleteSVD] = setup(matrix, svdArgs);

% Declare the intial structure
s = struct();

% Run the initial EOF on the Data
[s.eigVals, s.modes, s.Datax0, s.C] = simpleEOF(Data, matrix, svdArgs{:});

% Get the signals
s.signals = getSignals(s.Datax0, s.modes);

% Get the scaled signals
s.scaledSignals = scaleSignals(s.signals, s.eigVals);

% Get the explained variance
if incompleteSVD   % Use the data variance if svd was incomplete
    s.expVar = explainedVar(s.eigVals, s.Datax0);
else
    s.expVar = explainedVar(s.eigVals);
end    

% Rule N, and rotation require significance testing to continue
if ~blockMC

    % Run the Rule N significance test
    [s.nSig, s.randEigvals, s.normVals, s.thresh, s.trueConf, s.iterSigEigs, s.iterTrueConf] = ...
            ruleN(Data, matrix, s.eigVals, MC, noiseType, pval, svdArgs, showProgress, convergeTest);
    
    % Less than 2 significant eigs, cannot rotate
    if s.numSig < 2   
         % Do nothing
    
    % There are significant modes, rotate them...
    else 

        % Scale the eigenvectors
        s.scaModes = scaleModes( s.modes(:,1:s.nSig), s.eigVals(1:s.nSig));

        % Rotate the eigenvectors
        [s.rotModes, s.rotEigvals, s.rotExpVar, s.rotMatrix] = ...
            varimaxRotation( s.scaModes, s.modes(1:s.nSig));

        % Get the rotated signals
        s.rotSignals = getSignals(s.Datax0, s.rotModes);

        % Get the scaled rotated signals
        s.scaRotSignals = scaleSignals( s.rotSignals, s.rotEigvals);
    
    end

end
end

%%%%% Helper Functions %%%%%
function[showProgress, svdArgs, blockMC, convergeTest] = parseInputs(varargin)
inArgs = varargin;

% Set defaults
showProgress = 'noProgress';
svdArgs = {'svd'};
blockMC = false;
convergeTest = 'testConverge';

% Get input values
if ~isempty(inArgs)
    
    % Get each input
    for k = 1:length(inArgs)
        arg = inArgs{k};
        isSvdArg = false;
        
        if isSvdArg
            % Do nothing
        elseif strcmpi(arg, 'showProgress') 
            showProgress = 'showProgress';
        elseif strcmpi(arg, 'noSigTest')
            blockMC = true;
        elseif strcmpi(arg, 'noConvergeTest')
            convergeTest = 'noConvergeTest';
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
            
function[incompleteSvd] = setup(matrix, svdArgs)

incompleteSvd = false;

if length(svdArgs)
% Determine if explainedVar will need extra inputs
for ks = 1:length(varargin)
    spec = varargin{ks};
    
    if isscalar(spec)  % Only the first several eigs are stored
        incompleteSvd = true; 
    elseif strcmp(spec, 'econ') % Economy size SVD is performed
        incompleteSvd = true;
    end
    
    if notFullSVD && strcmp(covcorr, 'none')
        error('Explained variance cannot be calculated when an incomplete SVD is performed directly on a data matrix.');
    end
end
    
% SVDS is always economy sized
if strcmp(covcorr, 'svds')
    incompleteSvd = true;
end
    


end

