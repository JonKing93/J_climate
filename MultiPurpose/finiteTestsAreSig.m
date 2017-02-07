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
% nPass: The number of significance tests that pass the significance level p
%
% p: The signigicance levels. (e.g. p = 0.05 is 95% significant)
%
% ----- Outputs -----
%
% isSig: a boolean determining whether the tests still pass the
% significance level p
%
% nNeeded: The number of passed tests required to exceed the significance
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
errCheck(nTests, nPass, p);

% Calculate the binomial distribution probability for passing this many
% tests at the significance level.
prob = 0;
% For each test...
for k = nTests :-1: 0 
    % Get the probability of the number of passed tests (pnpt)    
    % (this is each term of the binomial expansion...)
    
    % Use exp(gammaln...) instead of nchoosek to avoid numerical precision errors
    % for large coefficients.
    pnpt = exp( gammaln(nTests+1) - gammaln(k+1) - gammaln(nTests-k+1)) * ...
        p^k * (1-p)^(nTests-k);
    
    % Get the sum to calculate the total probability.
    prob = pnpt + prob;
    
    % When the binomial probability exceeds the desired confidence level...
    if prob > p
        % Exit the loop
        break;
    end
end

% Move k to the previous value. This is number of tests expected to pass
% randomly.
nNeeded = k+1;

% Check the number of passed tests against k
% distribution
if nPass < nNeeded
    % The number of passed tests is less than the number randomly expected
    % to pass a finite number of tests with probability p. 
    
    % This dataset is not statistically significant at p
    isSig = false;
    
else
    % This many tests would not pass randomly. The tests are still
    % significant.
    isSig = true;
end

end

% ----- Helper Functions -----
function errCheck(nTests, nPass, p)

% Ensure p is between 0 and 1
if p < 0 || p > 1
    error('p must be between 0 and 1');
end

% Ensure test numbers are positive
if nTests < 1
    error('The number of tests must be positive');
end

if nPass < 0
    error('The number of passed tests cannot be negative');
end

% Ensure nPass is within the bounds of nTests
if nPass > nTests
    error('The number of passed tests cannot exceed the number of tests');
end

end