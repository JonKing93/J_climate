function[sigP, passTest] = fdr(p, q, varargin)
%% Controls the rate of false discovery for a set of p-values used to test
% multiple inter-related null hypotheses.
%
% [sigP, passTest] = fdr(p, q)
% Returns the (sorted) p values that are significant under the false
% discovery rate criteria for tests with any variable dependency.
%
% [sigP, passTest] = fdr(p, q, '+dep')
% Returns the (sorted) p values that are significant under the false
% discovery rate criteria for tests with ONLY positive dependency (i.e. all
% variables are positively correlated). This is a less conservative test
% than the case of any dependency.
%
%
% ----- Inputs -----
% 
% p: A set of p values from multiple null hypothesis tests. p may be
%       N-dimensional. NaN values are not given consideration in false 
%       discovery rate tests.
%
% q: The rate of false discovery. This is the percent of null hypotheses
%       that are falsely rejected and attributed as significant. q must be 
%       on the interval (0,1)
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
% A review on multiple comparisons including FDR:
%   Groppe, D.M., Urbach, T.P., & Kutas, M. (2011) Mass univariate analysis
%       of event-related brain potentials/fields I: A critical tutorial review. 
%       Psychophysiology, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x 


% Parse the inputs and do error checking
[posDependent] = parseInputs(varargin{:});
errCheck(p,q);

% Reshape p into a column vector
 pvec = p(:);

% Sort the p values
pvec = sort(pvec);

% Remove any NaN entries
pvec = pvec(~isnan(pvec));

% Get the total number of p values
m = numel(pvec);

% Assign an index i for each sorted p value
i = (1:m)';

% If the tests were on variables with ONLY positive dependency...
if posDependent
    
    % ... use the less conservative test: p <= i*q/m      
    passTest =  (pvec <= i*q/m);
    
else % Otherwise, use the more conservative test for any dependency
    denom = m * sum( 1./i );
    passTest =  (  pvec <=  i*q/denom  );
end
    

% Get the index of the last passed test
k = find(passTest, 1, 'last');

% Get the vector of significant p values
sigP = pvec(1:k);

% Get the boolean mask for the original p array
passTest =  (  p <= sigP(end)  );

end

% Helper Function
function[posDependent] = parseInputs( varargin )
inArgs = varargin;

% Set the default
posDependent = false;

if ~isempty(inArgs)
    if length(inArgs)==1 && strcmpi(inArgs{1}, '+dep')
        posDependent = true;
    else
        error('Unrecognized Input');
    end
end
end

function[] = errCheck(p, q)
if any( p(~isnan(p))>=1 | p(~isnan(p))<=0)
    error('All non-NaN p-values must be on the interval (0,1)');
elseif q<=0 || q>=1 || isnan(q)
    error('q must be on the interval (0,1)');
end
end
        
        
        
        
        
        
        
        
        
        
    
