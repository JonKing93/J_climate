function[corrmaps, pmaps, N] = lagFieldcorr( ts, iTS, field, iField, varargin)
%% Lags or leads a time series relative to a field and computes correlation
% coefficients and corresponding p-values.
% 
% [corrmaps, pmaps] = lagFieldcorr(ts, iTS, field, iField, 'lagTS', tsLags)
% Computes timeseries-field correlation coefficients and p-values at specifed
% lags and leads. Leads / lags the time series relative to the field.
%
% [corrmaps, pmaps] = lagFieldcorr(ts, iTS, field, iField, 'lagF', fLags)
% Computes timeseries-field correlation coefficients and p-values at specifed
% lags and leads. Leads / lags the field relative to the time series.
%
% [...] = lagFieldcorr(..., 'lagTS', tsLags, 'lagF', fLags)
% performs separate leads / lags for both the time series and the field. 
% 
% [...] = lagFieldcorr(..., 'fixedN')
% Restricts the size of the analyzed vector to ensure that correlations and 
% p-values are calculated using the same sample size for all leads and lags.
% Extra indices are removed from the end of the time series.
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
% tsLags: The desired lags (negative integers) and leads (positive integers)
%   at which to compute correlation coefficients. (e.g. lags = [-3 0 10]
%   will calculate the correlation with a lag of 3 time indices, no lags, and
%   lead of 10 indices for the time series.)
%
% fLags: The desired lags / leads for the field.
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
%   stored along the last dimension. Time series lags are recorded before
%   field lags.
%
% pmap: Displays the p value of the correlation coefficient for each point
%   on the field. Lags and Leads are stored along the last dimension
%
%
% ----- Written By -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu

% Parse inputs, error checking
[tsLags, fLags, fixN, iBounds, dim, corrArgs] = parseInputs(varargin{:});

% errCheck(iTS, iField, lags, iBounds);

% Get the indices of the points to use in each correlation
[tsPoints, fPoints] = getLagPoints(iTS, iField, tsLags, fLags, iBounds, fixN);

% Reshape the field into a 2D matrix along the dimension of interest
[field, dSize, dOrder] = dimNTodim2(field, dim);

% Preallocate the correlation and p value maps
nlags = numel(tsLags) + numel(fLags);
corrmaps = NaN( [1, size(field,2), nlags] );
pmaps = NaN( [1 size(field,2), nlags] );

% For each lead / lag...
for k = 1:nlags
    % Perform a field correlation
    [corrmaps(:,:,k), pmaps(:,:,k)] = ...
        fieldcorr( ts(tsPoints(:,k)), field(fPoints(:,k)), 'corrArgs', corrArgs);
end

% Reshape the maps to the original dimensionality
corrmaps = dim2TodimN( corrmaps, [1 dSize(2:end)], dOrder );
pmaps = dim2TodimN( pmaps, [1 dSize(2:end)], dOrder );

end

% % ----- Helper Functions -----
function[tsPoints, fPoints] = getLagPoints(iTS, iField, tsLags, fLags, iBounds)
% Get the 0-lag time intersection points
if isempty(iBounds)
    [~, xTS, xF] = intersect(iTS, iField);  % All intersecting time indices
else
    xTS = find( iTS>=iBounds(1)  &  iTS<=iBounds(2) );
    xF  = find( iField>=iBounds(1)  &  iField<=iBounds(2) );
end

% Preallocate logical arrays of the time indices to be compared for each lead / lag
nlags = numel(tsLags) + numel(fLags);
lTS = length(ts);
lF = size(field, 1);

tsPoints = false(lTS, nlags);
fPoints = false(lF, nlags);
lagcol = 1;  % Iterate through columns regardless of ts / field lag structure

% For the time series lags...
if ~isempty(tsLags)
    for k = 1:numel(tsLags)
       
        % Apply the lead / lag to the time series
        tsLagdex = xTS + tsLags(k);
        
        % Remove any values that have fallen off the time series time index
        tsLagdex( tsLagdex<=0 | tsLagdex>lTS ) = [];
        
        % Adjust the field indices if the time series is shorter
        fLagdex = xF;
        lDiff = numel(xTS) - numel(tsLagdex);
        
        if tsLags(k) > 0  
            % If the time series was lead over the edge, remove points from the end of the field      
            if lDiff > 0
                fLagdex(end-lDiff+1:end) = [];
            end

        elseif tsLags(k) < 0 
            % If the time series is lagged over the edge, remove points from the end of the field       
            if lDiff > 0
                fLagdex(1:lDiff) = [];
            end
        end
        
        % Set the desired indices in the logical arrays
        tsPoints(tsLagdex, lagcol) = true;
        fPoints(fLagdex, lagcol) = true;
        lagcol = lagcol + 1;
    end
end

% For the field lags
if ~isempty(fLags)
    for k = 1:numel(fLags)
        % Apply the lead / lag to the time series
        fLagdex = xF + fLags(k);
        
        % Remove any values that have fallen off the time series time index
        fLagdex( fLagdex<=0 | fLagdex>lF ) = [];
        
        % Adjust the field indices if the time series is shorter
        tsLagdex = xTS;
        lDiff = numel(xF) - numel(fLagdex);
        
        if fLags(k) > 0  
            % If the time series was lead over the edge, remove points from the end of the field      
            if lDiff > 0
                tsLagdex(end-lDiff+1:end) = [];
            end

        elseif fLags(k) < 0 
            % If the time series is lagged over the edge, remove points from the end of the field       
            if lDiff > 0
                tsLagdex(1:lDiff) = [];
            end
        end
        
        % Set the desired indices in the logical arrays
        tsPoints(tsLagdex, lagcol) = true;
        fPoints(fLagdex, lagcol) = true;
        lagcol = lagcol + 1;
    end
end
end

function[tsLags, fLags, fixN, iBounds, dim, corrArgs] = parseInputs(varargin)
inArgs = varargin;

% Set defaults
tsLags = [];
fLags = [];
fixN = false;
iBounds = [];
dim = 1;
corrArgs = {};

if ~isempty(inArgs)
    % Flag switches
    istslag = false;
    isflag = false;
    isBound = false;
    isDim = false;
    isCorrArg = false;
        
    for k = 1:length(inArgs)
        arg = inArgs{k};

        if istslag
            tsLags = arg;
            istslag = false;
        elseif isflag
            fLags = arg;
            isflag = false;
        elseif isBound
            iBounds = arg;
            isBound = false;
        elseif isDim
            dim = arg;
            isDim = false;
        elseif isCorrArg
            corrArgs = arg;
            isCorrArg = false;
        
        elseif strcmpi(arg, 'lagTS')
            if length(inArgs) >= k+1
                istslag = true;
            else
                error('tsLags not specified'); 
            end
            
        elseif strcmpi(arg, 'lagF')
            if length(inArgs) >= k+1
                isflag = true;
            else
                error('fLags not specified'); 
            end
            
        elseif strcmpi(arg, 'fixedN')
            fixN = true;
        
        elseif strcmpi(arg, 'restrictBounds')
            if length(inArgs) >= k+1
                isBound = true;
            else
                error('iBounds not specified'); 
            end
            
        elseif strcmpi(arg, 'fieldDim')
            if length(inArgs) >= k+1
                isDim = true;
            else
                error('fieldDim not specified'); 
            end
            
        elseif strcmpi(arg, 'corrArgs')
            if length(inArgs) >= k+1
                isCorrArg = true;
            else
                error('corrArgs not specified'); 
            end
            
        else
            error('Unrecognized input');
        end
    end
end    
end
%     
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % Do some error checking, get some sizes
% function[iTS, iField] = setup(tsIndices, fieldIndices, lags, indexBounds, lagTS, restrict)
% 
% % Check that indices are correctly formatted and overlapping
% if ~isvector(tsIndices)
%     error('tsIndices must be a vector');
% elseif ~isvector(fieldIndices)
%     error('fieldIndices must be a vector');   
% elseif any(mod(tsIndices,1) ~= 0) || any( mod(tsIndices,1) ~= 0)
%     error('Indices must be integers');
% elseif any( mod(tsIndices(2:end)-tsIndices(1:end-1),1) ~= 0) || ...
%         any( mod(fieldIndices(2:end)-fieldIndices(1:end-1),1) ~= 0)
%     error('Indices must increment by 1');
% end
% 
% % Indices are formatted correctly, get the overlap
% [overlap, iTS, iField] = intersect(tsIndices, fieldIndices);
% if isempty(overlap)
%     error('tsIndices and fieldIndices have no overlapping values');
% end
% 
% % Restrict the overlap if indexBounds are specified
% if restrict    
%     % Check that the index bounds are acceptable
%     if  ~any( overlap == indexBounds(1) )
%         error('Lower indexBound is not within the overlapping time indices');
%     elseif ~any( overlap == indexBounds(2) )
%         error('Upper indexBound is not within the overlapping time indices');
%     end
%     % Limit the overlap
%     minDex = find( overlap == indexBounds(1) );
%     maxDex = find( overlap == indexBounds(2) ); 
%     overlap = overlap(minDex:maxDex);
%     
%     % Update the indices
%     [~,iTS] = intersect( tsIndices, overlap);
%     [~,iField] = intersect(fieldIndices, overlap);
% end
% 
% % Check that lags is a vector of integers
% if ~isvector(lags)
%     error('lags must be a vector');
% elseif any( mod(lags, 1) ~= 0)
%     error('lags must be integers');
% end
% 
% % Check that lags are not too large
% maxlag = 0;
% maxlead = 0;
% if min(lags) < 0
%     maxlag = min(lags);
% end
% if max(lags) > 0
%     maxlead = max(lags);
% end
% 
% if lagTS % For lag/leading the time series
%    if ~any(tsIndices == (min(overlap)+maxlag) )
%        error('Maximum lag exceeds the length of the time series');
%    elseif ~any( tsIndices == max(overlap)+maxlead)
%        error('Maximum lead exceeds the length of the time series.');
%    end
% else
%    if ~any(fieldIndices == (min(overlap)+maxlag) )
%        error('Maximum lag exceeds the length of the time series');
%    elseif ~any( fieldIndices == max(overlap)+maxlead)
%        error('Maximum lead exceeds the length of the time series.');
%    end
% end
% 
% 
%    
% end
