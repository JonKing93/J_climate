function[maxFreq, maxPeriod] = maxFreqPeriod(singVecs)
%% Obtains the frequency and associated period with maximum power from the raw periodogram
% of a time series
%
% [maxFreq, maxPeriod] = maxFreqPeriod(singVecs)
%
% ----- Inputs -----
% 
% singVecs: A set of singular vectors. Each column is one vector.
% 
% ----- Outputs -----
%
% maxFreq: The frequency with maximum power for each vector
%
% maxPeriod: The period with maximum power for each vector

% Set the frequency function
freqFunc = @periodogram;

% Error check and the number of vectors
[nvecs] = setup(singVecs);

% Determine lengths and get w
[nF1, w] = freqFunc(singVecs(:,1));
lnf = length(nF1); % Length of the normed freq

% Preallocate
normFreq = NaN(lnf, nvecs);
maxDex = NaN(nvecs,1);

% For each column vector
for k = 1:nvecs
    % Compute a periodogram of the vector
    [normFreq(:,k)] = freqFunc( singVecs(:,k) );
    
    % Get the maximum frequency
    maxima = find( normFreq(:,k) == max(normFreq(:,k)),1 );
    
    % Give a warning if multiple maxima occur
    if length( find(normFreq(:,k) == max(normFreq(:,k)))) > 1
        warning('Multiple Maxima in vector %i',k);
    end
    maxDex(k) = maxima;
end

% Convert normalized frequencies to cycles / (time = N)
freq = w' ./ (2*pi);
period = 1 ./ freq;

% Get the maximum periods and frequencies
maxPeriod = period(maxDex)';
maxFreq = freq(maxDex)';

end


% ----- Helper Functions -----
function[nvecs] = setup(ts)

if ~ismatrix(ts)
    error('X must be a 2D matrix of column vectors');
end

[nvecs] = size(ts,2);

end