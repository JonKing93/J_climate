function[varargout] = multCompSig( p, N, nPass )
%% Returns a boolean to determine whether a finite number of passed 
% significance tests are acutally significant at the given significance
% levels.
%
% [nNeeded] = multCompSig( p, N )
% Returns the number of passed tests needed to remain significant given the
% number of hypothesis tests.
%
% [areSig, nNeeded] = multCompSig( p, N, nPass ) 
% Returns whether a specific number of passed significance tests remains
% significant given the number of hypothesis tests.
%
%
% ----- Inputs -----
% 
% p: A vector containing the tested significance levels. (e.g. p = 0.05 is 95%
%   significant). Each p values must be on the interval (0,1)
%
% N: The number of hypothesis tests
%
% nPass: A vector containing the number of tests that pass the significance level
%
%
% ----- Outputs -----
%
% areSig: a boolean vector determining whether the tests still pass each
% significance level p
%
% nNeeded: The number of passed tests required to exceed each significance
% level.
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

% Parse inputs
if nargin == 3
    havePass = true;
else
    havePass = false;
    nPass = NaN;
end

% Initial error checking
errCheck(p, N, nPass);

nNeeded = NaN(length(p),1);
for i = 1:length(p)
    
    % Get each binomial coefficient and add to cumulative probability
    cumProb = 0;
    for k = N(i):-1:0
        
        % Get the binomial probability distribution for the number of passed tests
        pnpt = binopdf(k, N(i), p(i));
        
        % Add to total probability
        cumProb = pnpt + cumProb;
        
        % When the binomial probability exceeds the desired confidence level...
        if cumProb > p(i)
            % Exit the loop
            break;
        end
    end
        
    % Move k to the previous value, this is nNeeded for N and p
    nNeeded(i) = k+1;
end

% If the number of passed tests is given, test for significance
if havePass
    % Check if the number of passed tests remains significant
    isSig = (nPass >= nNeeded);
    varargout(1:2) = {isSig, nNeeded};
else
    varargout = {nNeeded};
end

end

% ----- Helper Functions -----
function[] = errCheck(p, N, nPass)

if ~isvector(p)
    error('p must be a vector');
elseif ~isvector(N)
    error('N must be a vector');
elseif length(p) ~= length(N)
    error('p and N must have the same length');
end 
if ~isnan(nPass)
    if ~isvector(nPass)
        error('nPass must be a vector');
    elseif length(p) ~= length(nPass)
        error('p and nPass must have the same length');
    elseif any(nPass>N)
        error('The number of passed tests cannot exceed the number of tests');
    elseif any(nPass<0 | mod(nPass,1)~=0)
        error('The number of passed tests must be a non-negative integer');
    end
end

% Ensure p is between 0 and 1
if any(p < 0 | p > 1)
    error('p values must be between 0 and 1');
elseif any( N<1 | mod(N,1)~=0 )
    error('N must be a positive integer');
end

end