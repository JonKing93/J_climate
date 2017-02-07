function[randTS] = randNoiseSeries(ts, nSeries, noiseType)
% Build a set of random time series with the appropriate noise type.

% Get the length of the time series
lts = length(ts);
    
switch noise
    
    case 'white'
        randTS = randn( lts, MC );
        
    case 'red'
        % Preallocate
        randTS = NaN( lts, MC);
        
        % Get the lag-1 autocorrelation
        ar1 = corr( ts(1:end-1), ts(2:end));
        
        % Initialize the first row
        randTS(1,:) = randn(1,MC);
        
        % Calculate autocorrelation through matrix. Add random noise
        for k = 1:lts-1
            randTS(k+1,:) = (ar1 * randTS(k,:)) + randn(1,MC);
        end
        
        % Standardize so later scaling is correct
        randTS = zscore(randTS);
        
    otherwise
        error('Unrecognized noise type');
end

end