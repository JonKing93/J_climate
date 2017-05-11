function[sigP, passTest] = fdr(pvals, q, testType)
%% Controls the rate of false discovery for a set of p-values used to test
% multiple null hypotheses.
%
% [sigP, passTest] = fdr(pvals, q, testType)
% Returns the (sorted) p values that are significant under the false
% discovery rate criteria for tests with any variable dependency using the
% desired false discovery rate test.
%
%
% ----- Inputs -----
% 
% pvals: A set of p values from multiple null hypothesis tests. p may be
%       N-dimensional. NaN values are not given consideration in false 
%       discovery rate tests.
%
% q: The rate of false discovery. This is the percent of null hypotheses
%       that are falsely rejected and attributed as significant. (q = 0.05
%       is a commonly used value.) q must be on the interval (0,1)
%
% testType: A flag for the desired fdr procedure
%       'BH': The Benjamini-Hochberg procedure. Guaranteed to control the 
%               false detection rate when all data are independent or
%               positively correlated. Moderately conservative.
%       'BY': The Benjamini-Yekutieli procedure. Guaranteed for any type of
%               data dependency. A very conservative test.
%       'BKY': Benjamini-Krieger-Yekutieli. Guaranteed only for independent
%               datasets. Uses a two-stage process to provide a less overly
%               conservative approach than BH and BY.
%
%
% ----- Outputs -----
%
% sigP: A sorted column vector of the p values that pass the false
%       discovery rate criterion.
%
% passTest: A boolean array with the same dimensions as p. True at
%       locations that pass the false discovery rate criterion.
%
%
% ----- References -----
%
% Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery rate:
%   A practical and powerful approach to multiple testing. Journal of the 
%   Royal Statistical Society, Series B (Methodological). 57(1), 289-300. 
%
% Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery 
%   rate in multiple testing under dependency. The Annals of Statistics. 
%   29(4), 1165-1188. 
%
% Benjamini, Y., Krieger, A. M., & Yekutieli, D. (2006). Adaptive linear 
%   step-up procedures that control the false discovery rate. Biometrika, 491-507.
%
% A review on multiple comparisons including FDR:
%   Groppe, D.M., Urbach, T.P., & Kutas, M. (2011) Mass univariate analysis
%       of event-related brain potentials/fields I: A critical tutorial review. 
%       Psychophysiology, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x 
%
%
% ----- Written By -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu

% Do error checking
errCheck(pvals,q);

% Reshape p into a column vector
pvec = pvals(:);

% Sort the p values
pvec = sort(pvec);

% Remove any NaN entries
pvec = pvec(~isnan(pvec));

% Get the total number of p values
m = numel(pvec);

% Assign an index i for each sorted p value
i = (1:m)';

% Perform the appropriate test...
if strcmpi(testType, 'BH')       % Benjamini-Hochberg
    passTest =  (pvec <= i*q/m);
    
elseif strcmpi(testType, 'BY')   % Benjamini-Yekutieli
    denom = m * sum( 1./i );
    passTest =  (pvec <= i*q/denom);
    
elseif strcmpi(testType, 'BKY')  % Benjamini-Krieger-Yekutieli
    % Get q'
    q1 =    q / (1+q);
    
    % Get an estimate of the number of false null hypotheses (p values that
    % pass the test)
    r1 =  sum( pvec <= i*q1/m );
    
    % If nothing passed, stop, no null hypotheses are rejected
    if r1 == 0
        passTest = false( size(pvec) );
        
    else % Otherwise continue to the second stage
    
        % Get q'' and run the test
        q2 = (m / (m - r1)) * q1;
        passTest =  (pvec <= i*q2/m);
    end
end

% Get the index of the last passed test
k = find(passTest, 1, 'last');

% Get the vector of significant p values
sigP = pvec(1:k);

% Get the boolean mask for the original p array
if isempty(sigP)
    passTest = false( size(pvals) );
else
    passTest =  (  pvals <= sigP(end)  );
end

end

% ----- Helper Function -----

function[] = errCheck(p, q)
if any( p(~isnan(p))>=1 | p(~isnan(p))<=0)
    error('All non-NaN p-values must be on the interval (0,1)');
elseif q<=0 || q>=1 || isnan(q)
    error('q must be on the interval (0,1)');
end
end
        
        
        
        
        
        
        
        
        
        
    
