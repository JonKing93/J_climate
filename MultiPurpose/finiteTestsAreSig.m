function[outArg] = finiteTestsAreSig( nTests, p, varargin)
%% Determines whether a finite number of passed significance tests remain
% significant at the tested significance level.
%
% [nNeeded] = finiteTestsAreSig(nTests, p)
% Returns the minimum number of passed tests required to maintain significance
% when performing nTests significance tests at the significance level p.
%
% [isSig] = finiteTestsAreSig(nTests, p, nPass)
% Returns a boolean for whether a given number of passed significance tests
% remain significant at the tested signifcance level.
%
%
% ----- Inputs -----
% 
% nTests: The number of significance tests performed at a significance level, p
%
% p: The significance level used in the tests. p must be on the interval (0,1).
%
% nPass: The number of tests that pass the significance level.
%
%
% ----- Outputs -----
%
% isSig: a boolean determining whether the tests, in total, still pass the
%       significance level p
%
% nNeeded: The minimum number of passed tests required to remain above the
%       significance level.
%
% ----- References -----
%
% Livezey, R. E., & Chen, W. Y. (1983). Statistical field significance and 
%   its determination by Monte Carlo techniques. Monthly Weather Review, 111(1), 46-59.

% Initial error checking, parsing inputs
[nPass] = setup(nTests, nPass, varargin{:});

% Get the cumulative Probability of the Number of Passed Tests
cumProb = 0;
for k = nTests :-1: 0 

    % Get each coefficient of the binomial expansion. Use exp(gammaln...) 
    % instead of nchoosek to avoid numerical precision errors for large coefficients.
    pnpt = exp( gammaln(nTests+1) - gammaln(k+1) - gammaln(nTests-k+1)) * ...
        p^k * (1-p)^(nTests-k);

    % Get the sum to calculate the total probability.
    cumProb = pnpt + cumProb;

    % When the binomial probability exceeds the desired significance level...
    if cumProb > p
        % Exit the loop
        break;
    end
end

% Move k to the previous value. This is the minimum number of passed tests 
% needed to maintain significance at the desired p value.
nNeeded = k+1;

% Return the desired output
if nargin == 2
    outArg = nNeeded;
elseif nargin == 3
    outArg = ( nPass >= nNeeded );
else
    error('Incorrect number of inputs');
end

end

% ----- Helper Functions -----
function[nPass] = setup(nTests, p, varargin)

% Ensure p is between 0 and 1
if p < 0 || p > 1
    error('p values must be between 0 and 1');
end

% Ensure test numbers are positive
if nTests < 1 || mod(nTests,1)~=0
    error('The number of tests must be a positive integer');
end

% If given, get nPass, otherwise set to NaN
if ~isempty(varargin)
    if length(varargin) > 1
        error('Too many inputs');
    else
        nPass = varargin{1};
    end
    
    % Error check nPass
    if any(nPass) < 0
        error('The number of passed tests cannot be negative');
    elseif mod(nPass,1)~=0
        error('The numbre of passed tests must be an integer');        
    end
else
    nPass = NaN;
end

end