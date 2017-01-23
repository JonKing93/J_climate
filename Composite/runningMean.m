function[ S, meanDates] = runningMean( ts, window, varargin)
%% Computes a running mean over an n-dimensional dataset with a moving window.
%
% [S, meanDates] = runningMean(data, window)
%
% [S, meanDates] = runningMean(data, window, windowType, dates)
% 
% [S, meanDates] = runningMean(..., DIM)
%
% [S, meanDates] = runningMean(..., NANFLAG)
%
%
% ----- Inputs -----
%
% ts: An n-dimensional dataset of numeric values.
%
% window: The length of the window in which the running mean is computed.
%         If no additional arguments are given, runnningMean uses window as
%         the number of elements in data over which to calculate the mean.
%
% *** Optional Inputs ***
%
% windowType: Used to use a window with particular time duration (as
%         opposed to number of elements. Use of the windowType switches
%         requires input of the date input. The mean dates will be returned
%         with function output.
%   'index': Slides a window over the specified number of in the data
%           This differs from the defualt settings because it requires a
%           set of dates.
%   'years': Slides a window of the desired number of years.
%   'months': The window counts months
%   'days': The window counts days.
%   'duration': The window slides over a user-defined duration
%
% dates: Required with use of the windowType switches. An array containing
%   the dates of observations in the time series. Dates must be either
%   in datenum, or datevec format. 
%
% DIM: A scalar, specifies the dimension over which to compute the running
%   means. (Default: 1)
%
% NANFLAG: a switch specifying how to treat NaN
%   'includenan': Includes NaNs in running means. All means computed over a
%       NaN value will return NaN.
%   'omitnan': Running means will only be computed over non-NaN elements.
%
% 
% ----- Outputs -----
% 
% S: The new, running mean, time series
%
% meanDates: Requires input of the "dates" argument. A vector of containing
%   the dates on which each running mean is centered. The first column
%   contains the dates as datenum. The remaining columns contain the date
%   in a datevec format.
%
%
% ----- Written By -----
%
% Jonathan King, January 2017
%

% Process inputs, do some error checking
[dim, nanflag] = setup( varargin );

% Permute the data such that the dimension of the running mean is first.
dimOrder = 1:ndims(ts);
dimOrder(dim) = 1;
dimOrder(1) = dim;
ts = permute( ts, dimOrder);

% Reshape the data into a 2D matrix
sData = size(Data);
ts = reshape(ts, sData(1), prod(sData(2:end)));

% Get the indices over which the window will sit. (The window may not be
% fixed length if the data does not have an equal timestep -- This is often
% the case for monthly recorded data etc.)
switch windowType
    
    % This is the easy case, just move the window over the time series
    case index
        % Preallocate
        S = NaN( sData(1)-window+1, prod(sData(2:end)) );
        meanDates = NaN( sData(1)-window+1 );
        
        % Take the running mean
        for j = 1: (sData(1)-window+1)
            S(j) = mean( ts(j:j+window-1,:), nanflag );
            meanDates = mean( dates(j:j+window-1), nanflag );
        end
        
    case duration
        
        % Get a set of times with the added duration














end

% ----- Helper Functions -----
function[dates, dim, nanflag] = setup(inargs)

% Set default args
nanflag = 'includenan';
dim = 1;
dates = [];

% Set a few switches
dimSet = false;
isDates = false;
windowSet = false;
flagset = false;

% Deal with the optional input arguments
if ~isempty(inargs)
    
    % Get the number of input args
    nargs = length(inargs);
    
    % Process each input
    for k = 1:nargs
        
        % handles the windowType arguments
        if ~windowSet && isstring( inargs{k}  )
            switch inargs{k}
                case 'index'
                    
                case 'years'
                    
                case 'months'
                    
                case 'days'
                    
                case 'duration'
                    
                otherwise
                    if flagSet
                        error('Unrecognized input string');
                    end
            end
            isDates = true;
            windowSet = true;
            
        % handles the input dates
        elseif isDates
            dates = inargs{k};
            isDates = false;
        
        % Set the dimension on which to compute means
        elseif ~dimSet && isscalar(inargs{k})
            dim = inargs{k};
            dimSet = true;
            
        % set the nanflag
        elseif ~flagset && isstring( inargs{k} )
            switch inargs{k}
                case 'includenan'
                    nanflag = inargs{k};
                    
                case 'omitnan'
                    nanflag = inargs{k};
                    
                otherwise
                    error('Unrecognized input string');
            end
            flagset = true;
            
        else
            error('Unrecognized input');
        end
        
    end
            
end 
        
end
  