function[maxFreq, maxPeriod] = maxFreqPeriod(X)
%% Obtains the frequency and associated period with maximum power from the periodogram
% for each column of a set of column vectors
%
% -------------------------------------------------------------------------
% INPUTS:
% 
% X: a set of column vectors, must be 2D
% 
% -------------------------------------------------------------------------
% OUTPUTS:
%
% maxFreq: The frequency with maximum power for each vector
%
% maxPeriod: The period with maximum power for each vector

% Set the frequency function
freqFunc = @pmtm;


% Error check and get sizes
[nvecs] = setup(X);

% Determine lengths and get w
[nF1, w] = freqFunc(X(:,1));
lnf = length(nF1); % Length of the normed freq

% Preallocate
normFreq = NaN(lnf, nvecs);
maxDex = NaN(nvecs,1);

% Run over every column vector
for k = 1:nvecs
    % Compute a periodogram
    [normFreq(:,k)] = freqFunc( X(:,k) );
    
    % Get the maximum frequency
    maxima = find( normFreq(:,k) == max(normFreq(:,k)) );
    
    % Check for multiple maxima
    if length(maxima) > 1
        error('Multiple Maxima in vector %i',k);
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

function[nvecs] = setup(X)

if ~ismatrix(X)
    error('X must be a 2D matrix of column vectors');
end

[nvecs] = size(X,2);

end