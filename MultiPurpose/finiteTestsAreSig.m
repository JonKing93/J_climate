function[isSig, nNeeded] = finiteTestsAreSig( nTests, nPass, p )
%% Returns a boolean to determine whether a finite number of passed 
% significance tests are acutally significant at the given significance
% levels.
%
% [isSig] = finiteTestsAreSig( nTests, nPass, p )
%
% ----- Inputs -----
% 
% nTests: The number of significance tests
%
% nPass: A vector conatining the number of significance tests that pass 
%   each significance level p
%
% p: A vector containing the significance levels. (e.g. p = 0.05 is 95%
%   significant). Each p values must be on the interval (0,1)
%
% ----- Outputs -----
%
% isSig: a boolean vector determining whether the tests still pass each
% significance level p
%
% nNeeded: The number of passed tests required to exceed each significance
% level.
%
% ----- Sources -----
%
% Livezney and Chen (1983)
%
% ----- Author -----
%
% Jonathan King, 2017

% Initial error checking
[npvals] = setup(nTests, nPass, p);

% Calculate the binomial distribution probability for passing this many
% tests at the significance level.
pcount = 1;
nNeeded = NaN(npvals,1);

% For each pvalue
for q = p

    % Get the probability for each number of passed tests (pnpt). Each
    % value pnpt is one term in the binomal probability expansion.
    prob = 0;
    for k = nTests :-1: 0 
        % Get the probability of the number of passed tests (pnpt)    
        % (this is each term of the binomial expansion...)

        % Use exp(gammaln...) instead of nchoosek to avoid numerical precision errors
        % for large coefficients.
        pnpt = exp( gammaln(nTests+1) - gammaln(k+1) - gammaln(nTests-k+1)) * ...
            q^k * (1-q)^(nTests-k);

        % Get the sum to calculate the total probability.
        prob = pnpt + prob;

        % When the binomial probability exceeds the desired confidence level...
        if prob > q
            % Exit the loop
            break;
        end
    end

    % Move k to the previous value. This is the number of passed tests 
    % needed to maintain significance at the desired p value.
    nNeeded(pcount) = k+1;
    pcount = pcount + 1;
    
end

% Check if the number of passed tests remains significant
isSig = (nPass >= nNeeded);

end

% ----- Helper Functions -----
function[npval] = setup(nTests, nPass, p)

% Check p is a vector
if ~isvector(p)
    error('p must be a vector');
end

% Ensure p is between 0 and 1
if any(p < 0 | p > 1)
    error('p values must be between 0 and 1');
end

% Ensure test numbers are positive
if nTests < 1
    error('The number of tests must be positive');
end

% Check that nPass is a vector
if ~isvector(nPass)
    error('nPass must be a vector');
end

if any(nPass) < 0
    error('The number of passed tests cannot be negative');
end

% Ensure nPass is within the bounds of nTests
if any(nPass) > nTests
    error('The number of passed tests cannot exceed the number of tests');
end

% Check that nPass and p are the same length
lnPass = length(nPass);
npval = length(p);
if lnPass ~= npval
    error('The nPass and p vectors must be the same length');
end

end