function[randTS] = randNoiseSeries(ts, nSeries, noiseType)
%% Build a matrix of random time series with desired noise properties given
% an initial time series. All randomly generated time series will have a 
% of 0 and variance of 1.
%
% [randTS] = randNoiseSeries(ts, nSeries, noiseType)
%
% ----- Inputs -----
%
% ts: The time series on which to base the random time series. 
%
% nSeries: The number of random time series to generate
%
% noiseType: The type of noise to introduce into the random time series
%   'white': white Gaussian noise
%   'red': lag-1 correlated red noise
%
%
% ----- Output -----
%
% randTS: The matrix of randomly generated, noisy time series. Each column
%   of the matrix is a separate time series. Each time series has a mean of
%   0 and variance of 1.

% Do some quick error checking, make ts a column vector
ts = setup(ts, nSeries);

% Get the length of the time series
lts = length(ts);
    
switch noiseType
    
    case 'white'
        randTS = randn( lts, nSeries );
        
    case 'red'
        % Preallocate
        randTS = NaN( lts, nSeries);
        
        % Get the lag-1 autocorrelation
        ar1 = corr( ts(1:end-1), ts(2:end));
        
        % Initialize the first row
        randTS(1,:) = randn(1,nSeries);
        
        % Calculate autocorrelation through matrix. Add random noise
        for k = 1:lts-1
            randTS(k+1,:) = (ar1 * randTS(k,:)) + randn(1,nSeries);
        end
        
        % Standardize so later scaling is correct
        randTS = zscore(randTS);
        
    otherwise
        error('Unrecognized noise type');
end

end

% ----- Helper Functions -----
function[ts] = setup(ts, nSeries)

% Ensure ts is a column vector
if ~isvector(ts)
    error('ts must be a vector');
end
if isrow(ts)
    ts = ts';
end

% Ensure that nSeries is a scalar, positive integer
if ~isscalar(nSeries)
    error('nSeries must be a scalar, positive integer');
end
if mod(nSeries, 1) ~= 0
    error('nSeries must be an integer');
elseif nSeries < 1
    error('nSeries must be positive');
end
end
