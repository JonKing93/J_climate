function[areSig, nNeeded, nPassed] = pvalMultComp( p, pvals, d)
%% Checks a collection of p-values from mutliple hypothesis testing to check
% if they remain above a given significance level given test multiplicity
% 
% [areSig, nNeeded, nPassed] = pvalMultComp( p, pvals )
% Checks if a set of p values from multiple hypothesis tests remains above
% the desired significance levels.
% 
% [areSig, nNeeded, nPassed] = pvalMultComp( p, pvals, d )
% Checks if a multiple sets of p values from multiple hypothesis tests
% remain above the desired significance levels. Separate sets of p values
% are along dimension d of pvals.
%
%
% ----- Inputs -----
% 
% p: A vector of significance levels used for hypothesis testing. Elements 
%       of p must be on the interval (0,1). Each row of output will
%       correspond to a significance level.
%
% pvals: A set of p-values from multiple hypothesis tests. pvals may be
%       n-dimensional. If pvals contains multiple sets of p-values, the 
%       subsets must run along a common dimension. Elements of pvals must
%       be on the interval (0,1).
%
% d: Specifies a dimension of pvals that denotes separate sets of multiple
%       hypothesis tests. Each set of p-values will be individually tested
%       for significance. Each column of output will correspond to a
%       different set of p-values.
%
%
% ----- Outputs -----
%
% areSig: a logical specifying whether a set of multiple p-values remain 
%       above given significance level.
%
% nNeeded: The minimum number of passed tests required to remain above the
%       significance level.
%
% nPassed: The number of hypothesis tests initially passed by each set of
%       p-values.
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
if nargin == 2
    singleSet = true;
    d = NaN;
else
    singleSet = false;
end

% Error check, make everything a column vector
[p] = errCheck(p, pvals, d);

% Make the indexing dimension first for multiple sets of p-values
if ~singleSet
    dOrder = 1:ndims(pvals);
    dOrder(1) = d;
    dOrder(d) = 1;
    pvals = permute(pvals, dOrder);
end

% Get the number of hypothesis tests (N) for each set of p-values
if singleSet
    N = repmat( numel(pvals(~isnan(pvals))), size(p) );
    nSubsets = 1;
else
    % Get N for individual subsets
    nSubsets = size(pvals,1);
    N = NaN(length(p), nSubsets);
    
    for k = 1:nSubsets
        setK = pvals(k,:);
        N(:,k) = numel( setK(~isnan(setK)) );
    end
end

% Get the number of tests needed for each significance level and N
nNeeded = NaN(size(N));
for k = 1:nSubsets
    nNeeded(:,k) = multCompSig(p, N(:,k));
end

% Check how many tests are passed
nPassed = NaN( size(N) );
if singleSet
    for k = 1:length(p)
        pass =  ( pvals(~isnan(pvals)) <= p(k) );
        nPassed(k) = sum(pass(:));
    end
else
    for j = 1:nSubsets
        mapJ = pvals(j,:);
        for k = 1:length(p)
            passJ =   ( mapJ(~isnan(mapJ)) <= p(k) );
            nPassed(k,j) = sum(passJ(:));
        end
    end
end
    
% Check if the number of passed tests is sufficient to maintain
% significance
areSig =  (nPassed >= nNeeded);     

end               

function[p] = errCheck(p, pvals, d)
if any(p < 0 | p > 1)
    error('Significance levels (p) must be between 0 and 1');
elseif any(pvals(~isnan(pvals))<0 | pvals(~isnan(pvals))>1)
    error('pvals must be on the interval (0,1)');
end

if ~isnan(d) && (d>ndims(pvals) || d<0 || mod(d,1)~=0)
    error('d is not a valid dimension of pvals');
end

if isrow(p)
    p = p';
end
end