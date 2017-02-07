function[runMean] = runningMean(data, window, varargin)
%% This function gets the composite (discrete average) value of an n-dimensional 
% dataset along a dimension of choice.
%
%
% [runMean] = runningMean(data, window)
%
% [...] = compositeMean(..., dim)
%
% [...] = compositeMean(..., nanflag)
%
%
% ----- Inputs -----
%
% data: An n-dimensional dataset
%
% window: The number of indices over which to calculate each running mean.
%
% *** Optional Inputs ***
% 
% dim: (Default = 1) The dimension along which to take the mean
%
% nanflag: (Default = 'includenan') Specfies how to treat NaN values
%   'includenan': includes NaN values in mean
%   'omitnan': Do not include NaN values in mean
%
% ----- Outputs -----
%
% runMean: The running mean series.

% Initial error checking and sizing
[dim, nanflag, ldata] = setup(data, window, varargin);

% Convert data to 2D matrix with dimension of interest first.
[data, sData, dimOrder] = dimNTodim2(data, dim);

% Preallocate the composite mean
nMeans = ldata - window + 1;
runMean = NaN( nMeans, size(data,2) );

% Take the composite means
for k = 1:nMeans
    runMean(k,:) = mean( data(k:k+window-1,:), nanflag);
end

% Return to N-dimensional ordered dataset
sData(1) = nMeans;
runMean = dim2TodimN( runMean, sData, dimOrder);

end


% ----- Helper Functions -----
function[dim, nanflag, ldata] = setup(data, window, inArgs)

% Set defaults
dim = 1;
nanflag = 'includenan';

if length(inArgs) > 2
    error('Too many input arguments');
end
if ~isempty(inArgs)
    nanset = false;
    
    for k = 1:length(inArgs)
        if ~nanset && strcmpi(inArgs{k}, 'omitnan')
            nanset = true;
            nanflag = 'omitnan';
        elseif ~nanset && strcmpi(inArgs{k}, 'includenan')
            nanset = true;
            nanflag = 'includenan';
        elseif isscalar( inArgs{k} )
            dim = inArgs{k};
        else
            error('Unrecognized input');
        end
    end
    
end

% Check that the window is not longer than the data
ldata = size(data, dim);
if window > ldata
    error('Window length is longer than the data length');
end

end
        
        
        
        
        




