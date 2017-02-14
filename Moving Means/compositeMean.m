function[compMean] = compositeMean(data, window, varargin)
%% This function gets the composite (discrete average) mean of an n-dimensional 
% dataset along a dimension of choice.
%
% [compMean] = compositeMean(data, window)
% computes a composite mean along an N-dimensional dataset.
%
% [...] = compositeMean(..., skipObs)
% computes the composite mean and skips extra observations at either the
% start or end of the dataset
%
% [...] = compositeMean(..., dim)
% computes the composite mean along the a specified dimension.
%
% [...] = compositeMean(..., nanflag)
% computes the composite mean while treating NaNs in a user-specified
% manner
%
% ----- Inputs -----
%
% data: An n-dimensional dataset. Observations are assumed to be equally
%   spaced.
%
% window: The number of indices over which to calculate each composite
%   mean.
%
% *** Optional Inputs ***
% 
% skipObs: A flag specifying whether to skip extra observations from the
% beginning or end of the dataset. (By default, the number of observations
% in the dataset must be an exact multiple of window size. Specifying this
% flag loosens the restriction).
%       'skipStart': Skip extra observations at the beginning of the time series
%       'skipEnd': Skip extra observations at the end of the time series
%
% dim: (Default = 1) The dimension along which to take the mean
%
% nanflag: (Default = 'includenan') Specfies how to treat NaN values
%   'includenan': includes NaN values in mean
%   'omitnan': Do not include NaN values in mean
%
% ----- Outputs -----
%
% compMean: The composite mean series.

% Initial error checking and sizing
[dim, nanflag, ldata, skipObs] = setup(data, window, varargin);

% Convert data to 2D matrix with dimension of interest first.
[data, sData, dimOrder] = dimNTodim2(data, dim);

% Remove extra observations if required
extraObs = mod(ldata, window);
if extraObs>0
    if strcmpi(skipObs, 'begin')
        data = data( extraObs+1:end,:);
    elseif strcmpi(skipObs, 'end')
        data = data(1:end-extraObs,:);
    end
    ldata = size(data,1);
end

% Preallocate the composite mean
nComps = ldata / window;
compMean = NaN( nComps, size(data,2) );

% Take the composite means
for k = 1:nComps
    compMean(k,:) = mean(data(k*window-window+1:k*window,:), nanflag);
end

% Return to N-dimensional ordered dataset
sData(1) = nComps;
compMean = dim2TodimN( compMean, sData, dimOrder);

end


% ----- Helper Functions -----
function[dim, nanflag, ldata, skipObs] = setup(data, window, inArgs)

% Set defaults
dim = 1;
nanflag = 'includenan';
skipObs = NaN;

if length(inArgs) > 3
    error('Too many input arguments');
end
skipset = false;
if ~isempty(inArgs)
    nanset = false;
    dimset = false;
    
    for k = 1:length(inArgs)
        if ~nanset && strcmpi(inArgs{k}, 'omitnan')
            nanset = true;
            nanflag = 'omitnan';
        elseif ~nanset && strcmpi(inArgs{k}, 'includenan')
            nanset = true;
            nanflag = 'includenan';
        elseif isscalar( inArgs{k} ) && ~dimset
            dim = inArgs{k};
            dimset = true;
        elseif strcmpi(inArgs{k},'skipStart') && ~skipset
            skipObs = 'begin';
            skipset = true;
        elseif strcmpi(inArgs{k}, 'skipEnd') && ~skipset
            skipObs = 'end';
            skipset = true;
        else
            error('Unrecognized input');
        end
    end
    
end

% Check the window length is appropriate
ldata = size(data, dim);
if window > ldata
    error('Window length is longer than the data length');
end

% If not skipping observations, check that the data length is a multiple of
% window length
if ~skipset &&  (mod(ldata, window)~=0)
    error('To take a composite mean, the length of the data must be an exact multiple of the window length.');
end

end
