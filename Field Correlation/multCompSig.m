function[areSig, nNeeded] = multCompSig( p, varargin)
%% Accounts for the effects of multiple comparisons to determine whether
% multiple passed significance tests remain significant at the tested significance level.
%
% 
% [areSig, nNeeded] = multCompSig( p, pvals )
% Checks if a set of p values from multiple hypothesis tests remains above
% the desired significance levels.
% 
% [areSig, nNeeded] = multCompSig( p, pvals, d )
% Checks if a multiple sets of p values from multiple hypothesis tests
% remain above the desired significance levels. Separate sets of p values
% run along dimension d of pvals.
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
% p: The significance levels used in the tests; p must be on the interval 
%       (0,1). If p is a vector, each significance level will be tested and
%       each row of output will correspond to a different level.
%
% pvals: A set of p-values from multiple hypothesis tests. pvals may be
%       n-dimensional. p-values must be on the interval (0,1).
%
% d: Specifies a dimension of pvals that denotes separate sets of multiple
%       hypothesis tests. Each set of p-values will be individually tested
%       for significance. Each column of output will correspond to a
%       different set of p-values.
%
% N: The total number of hypothesis tests. If N is a vector, p must be a
%       vector of equal length.
%
% nPass: The number of passed significance tests. If nPass is a vector, N
%       and p must be of equal length.
%
%
% ----- Outputs -----
%
% areSig: a logical specifying whether a set of multiple p-values remain 
%       a given significance level.
%
% nNeeded: The minimum number of passed tests required to remain above the
%       significance level.
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

% Parse inputs, error check
[pvals, d, noPvals, N, nPass] = parseInputs(p, varargin{:});
errCheck(p, pvals, d, noPvals, N, nPass);

% If pvals are given, get the value(s) of N, the total number of tests
if ~isnan(pvals)    
    if isnan(d)
        N = repmat( numel(pvals(~isnan(pvals))), size(p) );
    else
        [N, p, pvals] = getSubsetNp(pvals, d, p);
    end
end

% Get the number of tests needed for each significance level, given each N
nNeeded = NaN( size(N) );
for i = 1:numel(nNeeded)

    % Initialize cumulative probability
    cumProb = 0;
    
    % Get each coefficient of the binomial expansion.
    for k = N(i): -1: 0
        
        % Use exp(gammaln...) instead of nchoosek to avoid numerical precision errors for large coefficients.
        pnpt = exp( gammaln( N(i)+1) - gammaln(k+1) - gammaln(N(i)-k+1)) * ...
            p(i)^k * (1-p(i))^(N(i)-k);

        % Get the sum to calculate the total probability.
        cumProb = pnpt + cumProb;

        % When the binomial probability exceeds the desired significance level...
        if cumProb > p
            % Exit the loop
            break;
        end

    end
    
    % Previous value of k is the minimum number of passed tests needed for
    % significance
    nNeeded(i) = k+1;
end

% Get the number of passed tests if p-values were given
if ~isnan(pvals)
    nPass = NaN( size(N) );  % Preallocate
    for k = 1:size(N,2)
        
        % Get the appropriate set of p-values
        if isnan(d)
            pSet = pvals;
        else
            pSet = pvals(k,:);
        end
        
        % Compare against each value of p
        for i = 1:size(p,1)
            nPass(i,k) =  (  pSet(~isnan(pSet)) <= p(i)  );
        end
    end
end
        
% Check if the number of passed tests exceeds the number of tests needed
if noPvals && nargin==2
    areSig = nNeeded;
    nNeeded = NaN;
else
    areSig = (nPass >= nNeeded);
end
end

% ----- Helper Functions -----
function[N, p, pvals] = getSubsetNp(pvals, d, p)

% Get the indexing dimension on the first dimension
dOrder = 1:ndims(pvals);
dOrder(1) = d;
dOrder(d) = 1;
pvals = permute(pvals, dOrder);

% Get the number of non-NaN entries in each subset              
sP = size(pvals);
nSubsets = sP(d);
N = NaN( 1, nSubsets );
p = repmat( p, [1, nSubsets]);

for k = 1:nSubsets
    setK = pvals(k,:);
    N(k) = numel( setK(~isnan(setK)) );
end

N = repmat( N, [length(p),1]);
end                    

function[pvals, d, noPvals, N, nPass] = parseInputs(p, varargin)
inArgs = varargin;

% Set defaults
pvals = NaN;
d = NaN;
N = NaN;
nPass = NaN;

if ~isempty(inArgs)
    
    arg1 = inArgs{1};
    if ~any(  arg1(~isnan(arg1))<=0  |  arg1(~isnan(arg1))>=1  )
        noPvals = false;
        pvals = arg1;
        
    elseif ~any( arg1(~isnan(arg1) < 1) ) &&  isvector(arg1) && length(arg1)==length(p)
        N = arg1;
        noPvals = true;
        
    else
        error('Unrecognized second input');
    end
    
    if length(inArgs) == 2
        if noPvals
            nPass = inArgs{2};
        else
            d = inArgs{2};
        end
    elseif length(inArgs) > 2
        error('Too many inputs');
    end       
else
    error('Insufficient inputs');
end

end

function[] = errCheck(p, pvals, d, noPvals, N, nPass)

if p < 0 || p > 1
    error('Significance levels (p) must be between 0 and 1');
end

if noPvals
    if ~isnan(nPass) 
        if(~isvector(nPass) || length(nPass)~=N)
            error('nPass must be a vector with the same length as p and N');
        elseif any(nPass<0 || mod(nPass,1)~=0)
            error('nPass must be a non-negative integer');
        end
    end
    if any(N<1 || mod(N,1)~=0)
        error('N must be a positive integer');
    end
    
else
    if any( pvals(~isnan(pvals))<0) || any(pvals(~isnan(pvals))>1)
        error('All pvals must be between 0 and 1');
    elseif d>ndims(pvals) || d<0 || mod(d,1)~=0
        error('d is not a valid dimension of pvals');
    end
end
end