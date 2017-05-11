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






% Get the number of tests needed for each significance level and N
nNeeded = NaN(size(N));
for k = 1:nSubsets
    nNeeded(:,k) = multCompSig(p, N(:,k));
end


    


end               

