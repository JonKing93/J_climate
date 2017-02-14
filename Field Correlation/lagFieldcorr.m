function[corrmaps, pmaps, isSig] = lagFieldcorr( ts, tsIndices, field, fieldIndices, lags, pvals, varargin)
%% Computes lagged or leading correlation coefficients for a time series and a field.
% Also tests the significance of the correlations.
% 
% [corrmaps, pmaps, isSig] = fieldcorr(ts, tsIndices, field, fieldIndices, lags, pvals, MC, noiseType)
% computes timeseries-field correlation coefficients at specified lags and 
% and tests whether correlation coefficients remain significant at a given
% level given a finite number of tests AND spatial reduction of degrees of
% freedom in a highly correlated field.
%
% [...] = fieldcorr(ts, tsIndices, field, fieldIndices, lags, pvals, 'noSpatial')
% suppresses significance testing with respect to spatial reduction of
% degrees of freedom.
% 
% [...] = fieldcorr(..., indexBounds)
% specifies a subset of indices over which to perform calculations.
%
% [...] = lagFieldcorr(..., lagType)
% specifies whether to lag/lead the time series or the field. (By default,
% lags and leads are applied to the time series.)
%
% [...] = fieldcorr(..., fieldDim)
% performs calculations along a specified dimension of the field. 
%
%
% ----- Inputs -----
% 
% ts: a time series vector.
%
% tsIndices: a set of time indices that overlap with fieldIndices
%
% field: A field. The time series at individual points in the field should
%   have the same length as ts. Field may be n-dimensional. By default,
%   fieldcorr assumes that the time series is along the first
%   dimension. (See the fieldDim argument if this is not the case)
%
% fieldIndices: A set of time indices for the field. These must overlap
%    with the tsIndices.
%
% lags: The desired lags (negative integers) and leads (positive integers)
%   at which to compute correlation coefficients. Lags and leads will shift
%   ts and/or field along the appropriate indices. (e.g. lags = [-3 0 10]
%   will calculate the correlation with a lag of 3 indices, no lags, and
%   lead of 10 indices.
%
% lagType: A switch specifying whether to lag/lead the time series or  field.
%       'ts': Lag/Lead the time series
%       'field': Lag/Lead the field
%
% pvals: A vector containing the p values of interest. (e.g. p = 0.05 
%   corresponds to the 95% confidence level. Correlation coefficients are
%   tested for significance at pval levels using Monte Carlo sampling.
%
% 'noSpatial': A flag to suppress significance testing with respect to
%   spatial loss of degrees of freedom.
%
% MC: The number of Monte Carlo simulations to use in spatial significance
%   testing.
%
% noisetype: The type of noise in the Monte Carlo time series.
%    'white': White Gaussian noise
%    'red': Red lag-1 autocorrelated noise with added Gaussian white noise
%
% indexBounds: A vector of two integers. Restricts the calculations to a
%   specific subset of indices. If unspecified, correlations are calculated
%   over ALL overlapping indices +/- leads and lags. 
%   (e.g. indexBounds = [48, 109] will restrict calculations to time
%   indices from 48-109 +/- lags and leads)
%
% fieldDim: (Default = 1) A scalar, specifies the dimension of the field on
%   which to perform the field correlation.
%
%
% ----- Outputs -----
%
% corrmaps: A N+1D array displaying the correlation coefficients at each
%   point of the field for each lag/lead. Different Lags and leads are
%   stored along the last dimension.
%
% pmap: Displays the p value of the correlation coefficient for each point
%   on the field. Lags and Leads are stored along the last dimension
%
% isSig: A logical matrix corresponding to each input p-value. A value of
%   true indicates a correlation is significant at the corresponding.
%   Each column contains the values for a particular lag or lead.
%
%
% ----- Additional Reading -----
%
% Livezney and Chen (1983)
%
%
% -----
% Jonathan King, 2017

%% Parse inputs, error checking, preallocation sizes, etc.

% Parse the inputs
[inArgs, fieldDim, restrict, indexBounds, lagTS] = parseInputs( varargin{:});

% Error check, get overlapping time indices
[iTS, iField] = setup(tsIndices, fieldIndices, lags, indexBounds, lagTS, restrict);

% Reshape the field into a 2D matrix along the dimension of interest
[field, dSize, dOrder] = dimNTodim2(field, fieldDim);

% Preallocate the final arrays
sField = size(field);
nPvals = length(pvals);
nLags = length(lags);

corrmaps = NaN( 1, sField(2), nLags );
pmaps = NaN( 1, sField(2), nLags  );
isSig = NaN( nPvals, nLags);


%% Get the initial (non-lagged/non-leading) indices of interest




%% Run a field correlation for each lag and lead

for k = 1:nLags
    
    % Lag the time series...
    if lagTS
        [corrmaps(1,:,k), pmaps(1,:,k), isSig(:,k)] = ...
            fieldcorr( ts( iTS+lags(k) ), field(iField,:), pvals, inArgs{:});
        
    % Or lag the field
    else
        [corrmaps(1,:,k), pmaps(1,:,k), isSig(:,k)] = ...
            fieldcorr( ts(iTS), field(iField+lags(k),:) , pvals, inArgs{:});
    end
    
end

% Reshape corrmaps and pmaps into the dimensionality of the original field
dOrder = [dOrder, length(dOrder)+1];
dSize = [1, dSize(2:end), nLags];

corrmaps = dim2TodimN( corrmaps, dSize, dOrder);
pmaps = dim2TodimN( pmaps, dSize, dOrder);


end

% ----- Helper Functions -----

% Reads variable inputs
function[inArgs, fieldDim, restrict, indexBounds, lagTS] = parseInputs(varargin)

% Set defaults
fieldDim = 1;
restrict = false;
indexBounds = NaN;
lagTS = true;

% Set other inputs
if isempty(varargin)
    error('Insufficient input arguments');
end

% Get the fieldcorr arguments separately
if strcmpi( varargin{1}, 'noSpatial')
    inArgs = varargin(1);
    iArg = 2;
else
    inArgs = varargin(2);
    iArg = 3;
end

% Set the lag field corr variables
lagSet = false;
dimSet = false;
boundSet = false;

while length(varargin) >= iArg 
    arg = varargin{iArg};
    
    if strcmpi(arg, 'ts') && ~lagSet
        lagTS = true;
        lagSet = true;
    elseif strcmpi(arg, 'field') && ~lagSet
        lagTS = false;
        lagSet = true;
    elseif isscalar(arg) && ~dimSet
        fieldDim = arg;
        dimSet = true;
    elseif isvector(arg) && length(arg)==2 && ~boundSet
        indexBounds = arg;
        restrict = true;
        boundSet = true;
    else
        error('Unrecognized input');
    end
    iArg = iArg+1;
end

end
    
% Do some error checking, get some sizes
function[iTS, iField] = setup(tsIndices, fieldIndices, lags, indexBounds, lagTS, restrict)

% Check that indices are correctly formatted and overlapping
if ~isvector(tsIndices)
    error('tsIndices must be a vector');
elseif ~isvector(fieldIndices)
    error('fieldIndices must be a vector');   
elseif any(mod(tsIndices,1) ~= 0) || any( mod(tsIndices,1) ~= 0)
    error('Indices must be integers');
elseif any( mod(tsIndices(2:end)-tsIndices(1:end-1),1) ~= 0) || ...
        any( mod(fieldIndices(2:end)-fieldIndices(1:end-1),1) ~= 0)
    error('Indices must increment by 1');
end

% Indices are formatted correctly, get the overlap
[overlap, iTS, iField] = intersect(tsIndices, fieldIndices);
if isempty(overlap)
    error('tsIndices and fieldIndices have no overlapping values');
end

% Restrict the overlap if indexBounds are specified
if restrict    
    % Check that the index bounds are acceptable
    if  ~any( overlap == indexBounds(1) )
        error('Lower indexBound is not within the overlapping time indices');
    elseif ~any( overlap == indexBounds(2) )
        error('Upper indexBound is not within the overlapping time indices');
    end
    % Limit the overlap
    minDex = find( overlap == indexBounds(1) );
    maxDex = find( overlap == indexBounds(2) ); 
    overlap = overlap(minDex:maxDex);
    
    % Update the indices
    [~,iTS] = intersect( tsIndices, overlap);
    [~,iField] = intersect(fieldIndices, overlap);
end

% Check that lags is a vector of integers
if ~isvector(lags)
    error('lags must be a vector');
elseif any( mod(lags, 1) ~= 0)
    error('lags must be integers');
end

% Check that lags are not too large
maxlag = 0;
maxlead = 0;
if min(lags) < 0
    maxlag = min(lags);
end
if max(lags) > 0
    maxlead = max(lags);
end

if lagTS % For lag/leading the time series
   if ~any(tsIndices == (min(overlap)+maxlag) )
       error('Maximum lag exceeds the length of the time series');
   elseif ~any( tsIndices == max(overlap)+maxlead)
       error('Maximum lead exceeds the length of the time series.');
   end
else
   if ~any(fieldIndices == (min(overlap)+maxlag) )
       error('Maximum lag exceeds the length of the time series');
   elseif ~any( fieldIndices == max(overlap)+maxlead)
       error('Maximum lead exceeds the length of the time series.');
   end
end


   
end
