function[varargout] = MC_fieldcorr(ts, field, MC, noiseType, p, varargin)
%% Performs a Monte Carlo test on a field correlation to determine if 
% correlations remain significant above a given level after accounting for
% spatial interdependence of field sites.
%
% [nNeeded, iterNPassed, iterTrueConf] = MC_fieldcorr(ts, field, MC, noiseType, p)
% Computes the number of passed significance tests required for the field
% correlation to remain above the significance levels, p, given observed spatial
% interdependence. Also tests the convergence on the Monte Carlo method.
%
% [areSig, nNeeded, iterNPassed, iterTrueConf] = MC_fieldcorr(ts, field, MC, noiseType, p, nPass)
% Returns a logical matrix for whether a given number of passed tests remains
% above each significance level p. Also tests Monte Carlo convergence.
%
% [...] = MC_fieldcorr(..., convergeFlag)
% May be used to block recording of Monte Carlo convergence.
%
% [...] = MC_fieldcorr(..., 'fieldDim', d)
% performs calculations along a specified dimension of the field. By
% default, correlations will be calculated along the first dimension.
%
% [...] = MC_fieldcorr(..., 'corrArgs', {corrArgs})
% Performs the test using alternative parameters for the "corr" function. 
% See the help section of "corr" for details.
%
%
% ----- Inputs -----
%
% ts: a time series. This is a single vector.
%
% field: A field. The time series at individual points in the field should
%       have the same length as ts. Field may be n-dimensional. By default,
%       fieldcorr assumes that the time series runs parallel to the first dimension.
%
% MC: The number of Monte Carlo iterations to perform
%
% noiseType: The type of noise to use during the test
%       'white': white Gaussian noise
%       'red': lag-1 autocorrelated noise with added Gaussian white noise
%
% p: A vector of significance levels to test
%
% nPass: A vector of the number of significance tests that initially pass a
%       significance level p. nPass be the same length as p.
%
% d: A scalar specifying the dimension of the field on which to perform the
%       correlation test.
%
% corrArgs: A cell containing alternative parameters for the "corr"
%       function, see the "corr" help section for details.
%
% convergeFlag: Specifies whether or not to record Monte Carlo convergence
%       'convergeTest': (Default) Records convergence data
%       'noConvergeTest': Does not record convergence data
%
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
%
%
% ----- References -----
%
% Livezey, R. E., & Chen, W. Y. (1983). Statistical field significance and 
%   its determination by Monte Carlo techniques. Monthly Weather Review, 111(1), 46-59.
%
%
% ----- Written By -----
%
% Jonathan King, 2017, University of Arizona, jonkin93@email.arizona.edu

% Parse inputs, do an error check
[nPass, corrArgs, convergeTest, dim] = parseInputs(varargin{:});
errCheck(ts, field, MC, p, nPass);

% Convert field to 2D along the dimension of interest
[field, ~, ~] = dimNTodim2(field, dim);

% Create a set of random time series with the desired noise properties
randTS = randNoiseSeries(noiseType, ts, MC);

% Get the significance of the correlation of each random series with each
% point on the field. The actual correlation coefficient is not needed here
[~, mcPval] = corr(randTS, field, corrArgs{:});

% Get the number of correlations that pass the significance level for each
% series.
mcNPass = NaN(MC, length(p));
for k = 1:length(p)
    mcNPass(:,k) = nansum(mcPval < p(k), 2);
end

% If testing Monte Carlo convergence...
if convergeTest
    
    % Preallocate the convergence data
    iterNPassed = NaN(MC,length(p));
    iterTrueConf = NaN(MC,length(p));
    
    % For each iteration...
    for k = 1:MC
        
        % ... sort the number of passed tests for the first k series
        iterPass = sort( mcNPass(1:k,:) );
        
        % Get the threshold index of the first significant element for the
        % current iteration.
        iterThresh = ceil( (1-p)*k );
        iterTrueConf(k,:) = iterThresh ./ k;
        
        % Get the number of required passed tests for the current iteration
        for j = 1:length(p)
            iterNPassed(k,j) = iterPass(iterThresh(j),j);
        end
    end
    
    % Get the number of passed tests needed to maintain significance
    nNeeded = iterNPassed(end,:);
    
else % ... not testing convergence
    % Sort the number of passed tests and get the threshold index
    mcNPass = sort(mcNPass);
    thresh = ceil( (1-p)*MC );
    
    % Get the number of passed tests at each significance level
    nNeeded = NaN(length(p),1);
    for k = 1:length(p)
        nNeeded(k) = mcNPass(thresh(k),k);
    end
end

% Return the desired output variable
if isnan(nPass)
    varargout(1) = {nNeeded};
else
    if ~iscolumn(nPass)
        nPass = nPass';
    end
    if ~iscolumn(nNeeded)
        nNeeded = nNeeded';
    end    
    varargout(1:2) = {(nPass > nNeeded), nNeeded}; % include the areSig output
end

if convergeTest
    varargout(end+1:end+2) = {iterNPassed, iterTrueConf};
else
    varargout(end+1:end+2) = {NaN, NaN};
end
    
end

% ----- Helper Methods -----
function[] = errCheck(ts, field2D, MC, p, nPass)

if ~isvector(ts)
    error('ts must be a vector');
elseif MC<1 || mod(MC,1)~=0
    error('MC must be a positive integer');
elseif any(p<=0 | p>=1)
    error('p must be on the interval (0,1)');
elseif length(ts) ~= size(field2D,1)
    error('The time series must run parallel to the first dimension of field2D');
end

if ~isnan(nPass)
    if ~isvector(nPass)
        error('nPass must be a vector')
    elseif ~isequal( length(nPass), length(p) )
        error('nPass must be the same length as p');
    end
end

end

function[nPassed, corrArgs, convergeTest, dim] = parseInputs(varargin)
inArgs = varargin;

% Set defaults
nPassed = NaN;
corrArgs = {};
convergeTest = true;
dim = 1;

if ~isempty(inArgs)
    
    isCorrArg = false;
    isDim = false;
    
    for k = 1:length(inArgs)
        arg = inArgs{k};
        
        % The nPassed arg
        if k==1 && isvector(arg)
            nPassed = arg;
        elseif isDim
            dim = arg;
            isDim = false;
        elseif isCorrArg
            corrArgs = arg;
            isCorrArg = false;
        elseif strcmpi(arg, 'fieldDim')
            isDim = true;
        elseif strcmpi(arg, 'corrArgs')
            isCorrArg = true;
        elseif strcmpi(arg, 'noConvergeTest')
            convergeTest = false;
        elseif strcmpi(arg, 'convergeTest')
            % Do nothing
        else
            error('Unrecognized Input');
        end
    end
end
end