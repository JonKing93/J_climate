function[corrmap, pmap, isSig] = fieldcorr(ts, field, pvals, varargin)
%% Determines the correlation coefficients of a time series with a field, 
% and tests for significance of the correlations.
% 
% [corrmap, pmap, isSig] = fieldcorr(ts, field, pvals, MC, noiseType)
% computes field-timeseries correlation coefficients and tests whether 
% correlation coefficients remain significant at a given level given a
% finite number of tests AND spatial reduction of degrees of freedom in a
% highly correlated field.
%
% [...] = fieldcorr(ts, field, pvals, 'noSpatial')
% suppresses significance testing with respect to spatial reduction of
% degrees of freedom
%
% [...] = fieldcorr(..., fieldDim)
% performs calculations along a specified dimension of the field.  
%
% ----- Inputs -----
% 
% ts: a time series. This is a single vector.
%
% field: A field. The time series at individual points in the field should
%   have the same length as ts. Field may be n-dimensional. By default,
%   fieldcorr assumes that the time series is along the first
%   dimension. (See the fieldDim argument if this is not the case)
%
% pvals: A vector containing the p values of interest. (e.g. p = 0.05 
%   corresponds to the 95% confidence level. Correlation coefficients are
%   tested for significance at pval levels using Monte Carlo sampling.
%
% 'noSpatial': A flag to suppress significance testing with respect to
%   spatial loss of degrees of freedom.
%
% noisetype: The type of noise in the Monte Carlo time series.
%    'white': White Gaussian noise
%    'red': Red lag-1 autocorrelated noise with added Gaussian white noise
%
% MC: The number of Monte Carlo simulations to use in spatial significance
%   testing.
%
% fieldDim: (Default = 1) A scalar, specifies the dimension of the field on
%   which to perform the field correlation.
%
% ----- Outputs -----
%
% corrmap: An array displaying the correlation coefficients at each point
% on the field.
%
% pmap: Displays the p value of the correlation coefficient for each point
% on the field.
%
% isSig: A boolean vector corresponding to each input p-value. A value of
% true indicates a correlation is significant at the corresponding p-value.
%
% ----- Additional Reading -----
%
% Livezney and Chen (1983)
%
% -----
% Jonathan King, 2017


%% Perform some initial error checking, get the number of dimensions. Make
% the time series a column vector.

[testSpatial, noiseType, MC, fieldDim] = parseInputs( varargin);

[ts] = setup(ts, field, pvals, fieldDim);


%% Perform the actual correlation

% Reshape the field into a 2D matrix along the dimension of interest
[field, dSize, dOrder] = dimNTodim2(field, fieldDim);

% Perform the correlation
[corrmap, pmap] = corr(ts, field);

% Reshape the correlation map to the shape of the original field
corrmap = dim2TodimN( corrmap, [1 dSize(2:end)], dOrder);

% p values of NaN are not actual tests, just missing data. Thus, remove
% them from statistical considerations.
ptests = pmap(~isnan(pmap));

% Now reshape the map of p-values
pmap = dim2TodimN( pmap, [1 dSize(2:end)], dOrder);


%% Perform statistical tests


%% Test 1 - Do the finite number of tests cause true significance to fall 
% below the desired significance level?

nTests = size(ptests, 2);
npvals = length(pvals);
nPass = NaN(npvals, 1);

% Get the number of correlations that initially pass each desired
% significance level.
for k = 1:npvals
    nPass(k) = sum( ptests <= pvals(k) );
end

% Check if finite tests are still significant at each significance level
isSig = finiteTestsAreSig( nTests, nPass, pvals);

% If none of the p values are significant, stop. There is no need to
% account for spatial reduction of degrees of freedom because none of them
% are significant anyways.
if ~any( isSig )
    return;
end


%% Test 2 - Check if spatial correlation reduces significance below 
% desired levels

% Only run if desired
if testSpatial
    
    % Begin by creating a set of random time series with the desired noise type.
    randTS = randNoiseSeries(ts, MC, noiseType);

    % Scale to the standard deviation of the time series
    randTS = randTS * std(ts);

    % Get the correlation significance of each Monte Carlo time series with
    % the data from each point in the field
    [~, mcPvals] = corr(randTS, field);

    % For each significance level
    mcPass = NaN( MC, npvals);
    for k = 1:npvals
        % For each random time series, get the number of points on the field
        % that pass the significance level.
        mcPass(:,k) = nansum( mcPvals < pvals(k), 2);

        % Sort the number of passed tests in ascending order
        mcPass(:,k) = sort( mcPass(:,k) );

        % Get the index of the p largest number of passed tests
        sigDex = ceil( pvals(k)*MC );
        mcNPass = mcPass(sigDex, k);

        % Compare this number of passed tests to the number passed by the time
        % series.
        if nPass(k) <= mcNPass
            % If the number of passed tests does not exceed the number from
            % the Monte Carlo simulation, then the null hypothesis is not
            % disproved for the desired significance level
            isSig(k) = false;
        end
    end
        
end

end


% ----- Helper Functions -----

% Reads the input arguments
function[testSpatial, noiseType, MC, fieldDim] = parseInputs( inArgs )

% Get the number of args
nIn = length(inArgs);

% Set default fieldDim
fieldDim = 1;

% See if the spatial testing is suppressed
if nIn == 0
    error('Insufficient input arguments');
end
if strcmpi(inArgs{1}, 'noSpatial')
    testSpatial = false;
    noiseType = NaN;
    MC = NaN;
    
    % Check for a fieldDim argument
    if nIn >= 2
        fieldDim = inArgs{2};
        if nIn > 2
            error('Too many input arguments');
        end 
    end
    
elseif nIn < 4
    testSpatial = true;
    MC = inArgs{1};
    noiseType = inArgs{2};
    
    % Check for a fieldDim argument
    if nIn == 3
        fieldDim = inArgs{3};
    end
    
else
    error('Too many input arguments');
end  
end


% Error checking and sizes, makes ts a column vector
function[ts] = setup(ts, field, pvals, fieldDim)

% Ensure fieldDim is scalar
if ~isscalar(fieldDim)
    error('fieldDim must be a scalar');
end

% Ensure fieldDim does not exceed the dimensions in field
sField = size(field);
if fieldDim > length(sField)
    error('fieldDim is greater than the number of dimensions in field');
end

% Ensure ts and field have the correct dimensions
if ~isvector(ts)
    error('ts must be a vector');
end

% Ensure the field is at least a vector
if isscalar(field)
    error('field cannot be a scalar');
end

% Ensure that ts and field have the same lengths
lTS = numel(ts);

if lTS ~= sField(fieldDim)
    error('ts must have the same number of observations as the dimension of interest in field');
end

% Make ts a column vector
if isrow(ts)
    ts = ts';
end

% Ensure that pvals is a vector of values only between 0 and 1
if ~isvector(pvals)
    error('pvals must be a vector');    
% Ensure the values are between 0 and 1
elseif any(pvals<=0) || any(pvals>=1)
    error('P values must be on the interval (0,1)');
end

end