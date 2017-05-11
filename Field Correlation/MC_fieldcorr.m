function[varargout] = MC_fieldcorr(ts, field2D, MC, noiseType, p, varargin)
%% Performs a Monte Carlo test on a field correlation to determine if 
% correlations remain significant above a given level after accounting for
% spatial interdependence of field sites.
%
%
% [nNeeded, iterNPassed, iterTrueConf] = MC_fieldcorr(ts, field2D, MC, noiseType, p)
% Computes the number of passed significance tests required for the field
% correlation to remain above the significance levels, p, given observed spatial
% interdependence. Also tests the convergence on the Monte Carlo method.
%
% [areSig, nNeeded, iterNPassed, iterTrueConf] = MC_fieldcorr(ts, field2D, MC, noiseType, p, nPassed)
% Returns a boolean for whether the given number of passed tests remains
% above the significance level p. Also tests Monte Carlo convergence.
%
% [...] = MC_fieldcorr(..., 'corrArgs', {corrArgs})
% Performs the test using alternative parameters for the "corr" function. 
% See the help section of "corr" for details.
%
% [...] = MC_fieldcorr(..., convergeFlag)
% May be used to block recording of Monte Carlo convergence.
%
%
% ----- Inputs -----
%
% ts: A time series vector. Cannot contain NaN.
%
% field2D: A field. The time series at individual points in the field should
%       have the same length as ts. Field must be 2 dimensional, and ts
%       must run parallel to the first dimension.
%
% MC: The number of Monte Carlo iterations to perform
%
% noiseType: The type of noise to use during the test
%       'white': white Gaussian noise
%       'red': lag-1 autocorrelated noise with added Gaussian white noise
%
% p: The significance level that the significance tests are performed against
%
% nPassed: The number of significance tests that passed the p significance level
%
% corrArgs: A cell containing alternative parameters for the "corr"
%       function, see the "corr" help section for details.
%
% convergeFlag: Specifies whether or not to record Monte Carlo convergence
%       'convergeTest': (Default) Records convergence data
%       'noConvergeTest': Does not record convergence data
%
% ----- Outputs -----
%
% nNeeded: The number of passed tests required to maintain significance at
%       the p significance level.
%
% areSig: A boolean for whether the given number of passed tests maintains
%       significance at the p significance level.
%
% iterNPassed: The number of passed significance tests that must be
%       exceeded to maintain significance for each Monte Carlo iteration.
%
% iterTrueConf: The true confidence level used for each iteration.


% Parse inputs, do an error check
[nPassed, corrArgs, convergeTest] = parseInputs(varargin{:});
errCheck(ts, field2D, MC, p);

% Create a set of random time series with the desired noise properties
randTS = randNoiseSeries(noiseType, ts, MC);

% Get the significance of the correlation of each random series with each
% point on the field. The actual correlation coefficient is not needed here
[~, mcPval] = corr(randTS, field2D, corrArgs{:});

% Get the number of correlations that pass the significance level for each
% series.
nPass = nansum(mcPval < p, 2);

% If testing Monte Carlo convergence...
if convergeTest
    
    % Preallocate the convergence data
    iterNPassed = NaN(MC,1);
    iterTrueConf = NaN(MC,1);
    
    % For each iteration...
    for k = 1:MC
        
        % ... sort the number of passed tests for the first k series
        iterPass = sort( nPass(1:k) );
        
        % Get the threshold index of the first significant element for the
        % current iteration.
        iterThresh = ceil( (1-p)*k );
        iterTrueConf(k) = iterThresh / k;
        
        % Get the number of required passed tests for the current iteration
        iterNPassed(k) = iterPass(iterThresh);
    end
    
    % Get the number of passed tests needed to maintain significance
    nNeeded = iterNPassed(end);
    
else % ... not testing convergence
    % Sort the number of passed tests and get the threshold index
    nPass = sort(nPass);
    thresh = ceil( (1-p)*MC );
    
    % Get the number of passed tests at the significance level
    nNeeded = nPass(thresh);
end

% Return the desired output variable
if isnan(nPassed)
    varargout(1) = {nNeeded};
else
    varargout(1:2) = {(nPassed > nNeeded), nNeeded}; % include the areSig output
end

if convergeTest
    varargout(end+1:end+2) = {iterNPassed, iterTrueConf};
end
    
end

% ----- Helper Methods -----
function[] = errCheck(ts, field2D, MC, p)

if ~isvector(ts)
    error('ts must be a vector');
elseif ~ismatrix(field2D)
    error('field2D must be a 2D matrix');
elseif MC<1 || mod(MC,1)~=0
    error('MC must be a positive integer');
elseif p<=0 || p>=1
    error('p must be on the interval (0,1)');
elseif length(ts) ~= size(field2D,1)
    error('The time series must run parallel to the first dimension of field2D');
end

end

function[nPassed, corrArgs, convergeTest] = parseInputs(varargin)
inArgs = varargin;

% Set defaults
nPassed = NaN;
corrArgs = {};
convergeTest = true;

if ~isempty(inArgs)
    
    isCorrArg = false;
    
    for k = 1:length(inArgs)
        arg = inArgs{k};
        
        % The nPassed arg
        if k==1 && isscalar(k)
            if mod(arg,1)==0 && arg>=0
                nPassed = arg;
            else
                error('nPassed must be a non-negative integer');
            end
        elseif isCorrArg
            corrArgs = arg;
        elseif strcmpi(arg, 'corrArgs')
            isCorrArg = true;
        elseif strcmpi(arg, 'noConvergeTest')
            convergeTest = false;
        elseif strcmpi(arg, 'convergeTest')
            % Do nothign
        else
            error('Unrecognized Input');
        end
    end
end
end