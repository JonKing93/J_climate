function[varargout] = multCompSig( p, N, nPass)
%% Returns a boolean to determine whether a finite number of passed 
% significance tests are acutally significant at the given significance
% levels.
%
% [nNeeded] = multCompSig( p, N )
% Returns the number of passed tests needed to remain significant given the
% number of hypothesis tests.[nN
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
% N: The number of hypothesis tests. If N is scalar, the same test size
%   will be applied to all significance levels. If N is a matrix, each row
%   must correspond with a significance level (p).
%
% nPass: A matrix containing the number of tests that pass each significance level.
%   For p, an (m x 1) vector, and N, a (z x n) matrix, nPass must be an (m x n)
%   matrix.
%
%
% ----- Outputs -----
%
% areSig: a logical matrix determining whether the number of passed tests 
%   remains above the significance level p.
%
% nNeeded: The number of passed tests required to maintain significance
%   levels. If p is (m x 1) and N is (z x n), then nNeeded will be an 
%   (m x n) matrix.
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

% Initial error checking, get whether N is scalar
[nScalar] = errCheck(p, N, nPass);

% Get the number of passes needed for each significance level and N
if havePass
    nNeeded = NaN( size(nPass) );
else
    nNeeded = NaN(length(p),size(N,2));
end
if nScalar
    N = repmat(N, [length(p),1]);
end

% For each significance level and test size
for i = 1:length(p)
    for j = 1:numel(N(i,:))
        
        % Get each binomial coefficient and add to cumulative probability
        cumProb = 0; 
        for k = N(j):-1:0

            % Get the binomial probability distribution for the number of passed tests
            pnpt = binopdf(k, N(j), p(i));

            % Add to total probability
            cumProb = pnpt + cumProb;

            % When the binomial probability exceeds the desired confidence level...
            if cumProb > p(i)
                % ...exit the loop
                break;
            end
        end

        % Move k to the previous value, this is nNeeded for N and p
        if nScalar
            nNeeded(i,:) = k+1;
        else
            nNeeded(i,j) = k+1;
        end
    end
end

% If the number of passed tests is given, test for significance
if havePass
    % Check if the number of passed tests remains significant
    areSig = (nPass >= nNeeded);
    varargout(1:2) = {areSig, nNeeded};
else
    varargout = {nNeeded};
end

end

% ----- Helper Functions -----
function[nScalar] = errCheck(p, N, nPass)

if ~isvector(p)
    error('p must be a vector');
end
if isscalar(N)
    nScalar = true;
else
    nScalar = false;
    if ~ismatrix(N)
        error('N must be a matrix');
    elseif size(N,1) ~= length(p)
        error('p and N must have the same length');
    end
end
    
if ~isnan(nPass)
    if ~ismatrix(nPass)
        error('nPass must be a matrix');
    elseif size(nPass,1)~=length(p)
        error('nPass must have a row for each significance level (p)');
    elseif ~nScalar && ~isequal( size(N), size(nPass))
        error('N and nPass must have the same size');
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