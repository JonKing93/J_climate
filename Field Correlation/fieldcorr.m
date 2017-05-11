function[corrmap, pmap] = fieldcorr(ts, field, varargin)
%% Determines the correlation coefficients of a time series with a field,
% as well as the significance of the correlation.
% 
% [corrmap, pmap] = fieldcorr(ts, field)
% computes field correlations for a time series. Returns correlation
% coefficients and significance p values.% highly correlated field.
%
% [...] = fieldcorr(..., 'fieldDim', d)
% performs calculations along a specified dimension of the field. By
% default, correlations will be calculated along the first dimension.
%
% [...] = fieldcorr(..., 'corrArgs', {corrArgs})
% Computes correlations using alternative parameters for the "corr"
% function. See the help section of "corr" for details.
%
%
% ----- Inputs -----
% 
% ts: a time series. This is a single vector.
%
% field: A field. The time series at individual points in the field should
%       have the same length as ts. Field may be n-dimensional. By default,
%       fieldcorr assumes that the time series runs parallel to the first dimension.
%
% d: A scalar specifying the dimension of the field on which to perform the
%       correlation.
%
% corrArgs: A cell containing alternative parameters for the "corr"
%       function, see the "corr" help section for details.
%
%
% ----- Outputs -----
%
% corrmap: An array displaying the correlation coefficients at each point
%       on the field.
%
% pmap: Displays the p value of the correlation coefficient for each point
%       on the field.
%
%
% ----- Written By -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu


% Parse the inputs, do some error checking
[fieldDim, corrArgs] = parseInputs( varargin{:});
[ts] = errCheck(ts, field, fieldDim); % Also makes ts a column vector

% Reshape the field into a 2D matrix along the dimension of interest
[field, dSize, dOrder] = dimNTodim2(field, fieldDim);

% Perform the correlation
[corrmap, pmap] = corr(ts, field, corrArgs{:});

% Reshape the correlation map to the shape of the original field
corrmap = dim2TodimN( corrmap, [1 dSize(2:end)], dOrder);

% Reshape the map of p-values
pmap = dim2TodimN( pmap, [1 dSize(2:end)], dOrder);

end


% ----- Helper Functions -----

% Reads the input arguments
function[fieldDim, corrArgs] = parseInputs( varargin )
inArgs = varargin;

% Set defaults
fieldDim = 1;
corrArgs = {'type', 'Pearson'};

if ~isempty(inArgs)
    iscorrarg = false;
    isdim = false;
    
    for k = 1:length(inArgs)
        arg = inArgs{k};
        
        if iscorrarg
            corrArgs = arg;
            iscorrarg = false;
        elseif isdim
            fieldDim = arg;
            isdim = false;
        elseif strcmpi(arg, 'fieldDim')
            if length(inArgs) >= k+1
                isdim = true;
            else
                error('fieldDim was not specified');
            end
        elseif strcmpi(arg, 'corrArgs')
            if length(inArgs) >= k+1
                iscorrarg = true;
            else
                error('corrArg was not specified');
            end
        else
            error('Unrecognized input');
        end
    end
end
end

% Error checking and sizes, makes ts a column vector
function[ts] = errCheck(ts, field, fieldDim)

% Ensure fieldDim is scalar
if ~isscalar(fieldDim)
    error('fieldDim must be a scalar');
end

% Ensure fieldDim does not exceed the dimensions in field
sField = size(field);
if fieldDim > length(sField)
    error('fieldDim is greater than the number of dimensions in field');
end

% Ensure ts and field have the correct dimensions
if ~isvector(ts)
    error('ts must be a vector');
end

% Ensure that ts and field have the same lengths
lTS = numel(ts);
if lTS ~= sField(fieldDim)
    error('ts must have the same number of observations as the dimension of interest in field');
end

% Make ts a column vector
if isrow(ts)
    ts = ts';
end
end