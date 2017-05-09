function[corrmaps, pmaps] = lagFieldcorr( ts, iTS, field, iField, lags, varargin)
%% Lags or leads a time series relative to a field and computes correlation
% coefficients and corresponding p-values.
% 
% [corrmaps, pmaps, isSig] = lagFieldcorr(ts, iTS, field, iField, lags)
% Computes timeseries-field correlation coefficients and p-values at specifed
% lags and leadsof the time series.
% 
% [...] = lagFieldcorr(..., 'fixedN')
% Restricts the size of the analyzed vector to ensure that correlations and 
% p-values are calculated using the same sample size for all leads and lags.
%
% [...] = lagFieldcorr(..., 'restrictBounds', iBounds)
% Specifies a specific subset of indices over which a 0-lag correlation
% would be computed. All lags and leads are applied to the restricted
% time series vector.
%
% [...] = lagFieldcorr(..., 'fieldDim', d)
% performs calculations along a specified dimension of the field. By
% default, correlations will be calculated along the first dimension.
%
% [...] = lagFieldcorr(..., 'corrArgs', {corrArgs})
% Computes correlations using alternative parameters for the "corr"
% function. See the help section of "corr" for details.
%
%
% ----- Inputs -----
% 
% ts: a time series vector.
%
% iTS: a set of time indices for time series. See the assignLagIndices function.
%
% field: A field. The time series at individual points in the field should
%   have the same length as ts. Field may be n-dimensional. By default,
%   fieldcorr assumes that the time series is along the first
%   dimension.
%
% iField: A set of time indices for the field. These must overlap
%    with the tsIndices. See the assignLagIndices function.
%
% lags: The desired lags (negative integers) and leads (positive integers)
%   at which to compute correlation coefficients. (e.g. lags = [-3 0 10]
%   will calculate the correlation with a lag of 3 indices, no lags, and
%   lead of 10 indices for the time series.) Note that lagging the 
%   timeseries is equivalent to leading the field.
%
% indexBounds: A vector containing the minimum and maximum indices over
%   which to define a 0-lagged correlation.
%
% fieldDim: (Default = 1) A scalar, specifies the dimension of the field on
%   which to perform the field correlation.
%
% corrArgs: A cell containing alternative parameters for the "corr"
%       function, see the "corr" help section for details.
%
%
% ----- Outputs -----
%
% corrmaps: An (n+1)D array displaying the correlation coefficients at each
%   point of the field for each lag/lead. Different Lags and leads are
%   stored along the last dimension.
%
% pmap: Displays the p value of the correlation coefficient for each point
%   on the field. Lags and Leads are stored along the last dimension
%
%
% ----- Written By -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu

% Parse inputs, error checking, preallocation sizes, etc.
[fixN, iBounds, dim, corrArgs] = parseInputs(varargin{:});

% Error check, get overlapping time indices
[iTS, iField] = setup(iTS, iField, lags, indexBounds, lagTS, restrict);

% Reshape the field into a 2D matrix along the dimension of interest
[field, dSize, dOrder] = dimNTodim2(field, fieldDim);

% Preallocate the final arrays
sField = size(field);
nPvals = length(pvals);
nLags = length(lags);

corrmaps = NaN( 1, sField(2), nLags );
pmaps = NaN( 1, sField(2), nLags  );
isSig = NaN( nPvals, nLags);






%% Run a field correlation for each lag and lead

for k = 1:nLags
    k
    
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

function[fixN, iBounds, dim, corrArgs] = parseInputs(varargin)
inArgs = varargin;

% Set defaults
fixN = false;
iBounds = [];
dim = 1;
corrArgs = {}


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
