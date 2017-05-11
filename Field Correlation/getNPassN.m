function[nPass, N] = getNPassN(p, pvals, d)
%% Gets the total number of hypothesis tests that initially pass a
% significance level, along with the number of tests.
%
% [nPass, N] = getNPassN(p, pvals)
% Gets the total number of tests and passed tests for a set of p-values at
% different significance levels.
%
% [nPass, N] = getNPassN(p, pvals, d)
% Gets the total number of tests and passed tests for multiple set of
% p-values. Separate sets of p values are along dimension d.
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
% nPassed: The number of hypothesis tests initially passed by each set of
%       p-values.
%
% N: The total number of hypothesis tests for each set of p-values.
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

% Check how many tests are passed
nPass = NaN( size(N) );
if singleSet
    for k = 1:length(p)
        pass =  ( pvals(~isnan(pvals)) <= p(k) );
        nPass(k) = sum(pass(:));
    end
else
    for j = 1:nSubsets
        mapJ = pvals(j,:);
        for k = 1:length(p)
            passJ =   ( mapJ(~isnan(mapJ)) <= p(k) );
            nPass(k,j) = sum(passJ(:));
        end
    end
end

end

% ----- Helper Function -----
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