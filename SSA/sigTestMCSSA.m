function[sigEigDex, highSurrTail, lowSurrTail] = sigTestMCSSA( p, dataEigvals, surrEigvals)
%% Performs a significance test of data SSA eigenvalues against MC_SSA surrogates eigenvalues
%
% [sigEigDex, highSurrTail, lowSurrTail] = sigTestMCSSA( p, dataEigvals, surrEigvals)
%
% ----- Inputs -----
%
% p: The confidence interval that the test should pass. This should be
% a number between 0 and 1 (e.g. 0.95 is the 95% confidence interval)
%
% dataEigval: A set of time series eigenvalues. Each column contains the
% eigenvalues of a particular time series
%
% surrEigvals: A set of surrogate eigenvalues generated using MC_SSA. dim1 is
% the Monte Carlo number, dim2 is each windowed eigenvalue, dim3 is the
% time series.
%
% ----- Outputs -----
%
% sigEigDex: A boolean (true/false) vector with the indices of the data 
% eigenvalues that pass the significance test. Each column corresponds to a
% particular time series.
%
% highSurrTail: A matrix of the surrogate eigenvalues on the upper tail of
% the significance test. Each column corresponds to a particular time
% series.
%
% lowSurrTail: A matrix of the surrogate eigenvalues on the lower tail of
% the significance test. Each column corresponds to a particular time
% series.
%

% Run an error check and get some initial sizes.
[MC] = setup(p, dataEigvals, surrEigvals);

% Get the indices of the confidence interval
maxRowDex = ceil( MC * (p + (1 - p)/2)  ) ;
minRowDex = floor( MC * (1-p)/2 );

% Get the vectors of the two tails of the confidence interval;
highSurrTail = squeeze( surrEigvals( maxRowDex, :, :) );
lowSurrTail = squeeze( surrEigvals( minRowDex, :, :) );

% For a single time series, the squeeze operation will convert highSurrTail
% and lowSurrTail to row vectors. This returns them to column vectors.
if size(highSurrTail,1) == 1
    highSurrTail = highSurrTail';
    lowSurrTail = lowSurrTail';
end

% Get the indices of the eigenvalues that passed the significance test
sigEigDex = dataEigvals > highSurrTail;

end


% ----- Helper functions -----
function[MC] = setup(pval, dataEigval, surrEigval)

% Ensure pval is between 0 and 1
if pval >= 1 || pval <= 0
    error('Confidence interval must be between 0 and 1');
end

% Ensure the inputs have the correct dimensionality
if ~ismatrix(dataEigval)
    error('dataEigval must be a 2D matrix');
end

if ndims(surrEigval) > 3
    error('surrEigval cannot have more than 3 dimensions');
end

% Ensure the dimensions of the inputs align
[dwindow, dnseries] = size(dataEigval);
[MC, swindow, snseries] = size(surrEigval);

if dwindow ~= swindow || dnseries ~= snseries
    error('Data and Surrogate eigenvalue matrices are not correctly aligned');
end

% Ensure the desired confidence interval is possible for the MC number
if floor( MC * (1-pval)/2 ) == 0
    minMC = 1/((1-pval)/2);
    error('The Monte Carlo number is too small for this confidence interval. A value of %i or larger is required.',minMC);
end
end