function[corrmaps, pmaps, N, tsPoints, fPoints] = lagFieldcorr( ts, iTS, field, iField, varargin)
%% Lags or leads a time series relative to a field and computes correlation
% coefficients and corresponding p-values.
% 
% [corrmaps, pmaps, N, tsPoints, fPoints] = lagFieldcorr(ts, iTS, field, iField, 'lagTS', tsLags)
% Computes timeseries-field correlation coefficients and p-values at specifed
% lags and leads. Leads / lags the time series relative to the field.
%
% [...] = lagFieldcorr(ts, iTS, field, iField, 'lagF', fLags)
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
% N: The number of points used in each correlation calculation
%
% tsPoints: Indices to the time series points used in each correlation.
%   Each column represents a specific lead or lag. Time series leads/lags are
%   listed before field lags/leads.
%
% fPoints: Indices to the field points used in each correlation calculation.
%
%
% ----- Written By -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu

% Parse inputs, error checking
[tsLags, fLags, fixN, iBounds, dim, corrArgs] = parseInputs(varargin{:});
[ts, iTS, iField] = errCheck(ts, iTS, iField, tsLags, fLags, iBounds); % Also converts row vectors to column

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
        fieldcorr( ts(tsPoints(:,k)), field(fPoints(:,k),:), 'corrArgs', corrArgs);
end

% Reshape the maps to the original dimensionality
corrmaps = dim2TodimN( corrmaps, [1 dSize(2:end) nlags], [dOrder max(dOrder)+1] );
pmaps = dim2TodimN( pmaps, [1 dSize(2:end) nlags], [dOrder max(dOrder)+1] );

% Get the sample size for each set of correlations
N = sum(tsPoints)';

end

% % ----- Helper Functions -----
function[tsPoints, fPoints] = getLagPoints(iTS, iField, tsLags, fLags, iBounds, fixN)
% Get the 0-lag time intersection points
if isempty(iBounds)
    [~, xTS, xF] = intersect(iTS, iField);  % All intersecting time indices
else
    xTS = find( iTS>=iBounds(1)  &  iTS<=iBounds(2) );
    xF  = find( iField>=iBounds(1)  &  iField<=iBounds(2) );
end

% Preallocate logical arrays of the time indices to be compared for each lead / lag
nlags = numel(tsLags) + numel(fLags);
lTS = length(iTS);
lF = length(iField);

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

if fixN
    % Limit all comparisons to the same number of points
    N = min( sum(tsPoints) );

    % For each column of lag points
    for k = 1:nlags
        % Get the first N intersection points
        saveTS = find( tsPoints(:,k), N );
        saveF = find( fPoints(:,k), N );

        % Remove all other intersection points
        tsPoints(:,k) = false;
        fPoints(:,k) = false;

        tsPoints(saveTS, k) = true;
        fPoints(saveF, k) = true;
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
corrArgs = {'type', 'Pearson'};

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
    
function[ts, iTS, iField] = errCheck(ts, iTS, iField, tsLags, fLags, iBounds)

% Check time indices
if ~isvector(ts) || ~isvector(iTS) || ~isvector(iField)
    error('ts, iTS, and iField must be vectors');
elseif any(mod(iTS,1)~=0) || any(mod(iField,1)~=0)
    error('iTS and iField must be integers');
elseif any( mod(iTS(2:end)-iTS(1:end-1),1) ~= 0 ) || ...
        any( mod(iField(2:end)-iField(1:end-1),1) ~= 0)
    error('Indices must increment by 1');
elseif isempty( intersect(iTS, iField) )
    error('iTS and iField have no overlapping values');
end

% Check the lags
if ~isempty(tsLags)
    if ~isvector(tsLags)
        error('tsLags must be a vector')
    elseif any( mod(tsLags,1) ~= 0 )
        error('tsLags must be integers');
    elseif max(abs(tsLags)) >= length(iTS)
        error('A time series lag exceeds the time series length');
    end
end
if ~isempty(fLags)
    if ~isvector(fLags)
        error('fieldLags must be a vector')
    elseif any( mod(fLags,1) ~= 0 )
        error('fieldLags must be integers');
    elseif max(abs(fLags)) >= length(iField)
        error('A field lag exceeds the field length');
    end
end
if isempty(tsLags) && isempty(fLags)
    error('No leads or lags have been specified');
end
    
% Column vectors
if ~iscolumn(ts)
    ts = ts';
end
if ~iscolumn(iTS)
    iTS = iTS';
end
if ~iscolumn(iField)
    iField = iField';
end

% Check the iBounds
if ~isempty(iBounds)
    overlap = intersect(iTS, iField);
    if ~any( overlap == iBounds(1) )
        error('Lower iBound is not within overlapping time indices');
    elseif ~any( overlap == iBounds(2) )
        error('Upper iBound is not within overlapping time indices');
    end
end  
end
