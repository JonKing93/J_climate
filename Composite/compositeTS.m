function[compTS, compDatenum, compDatevec] = compositeTS(ts, obsDatevec, obsDatenum, monthYear, nYears, varargin)
%% This function gets the composite (average) value of an n-dimensional 
% dataset along the first dimension.
%
% monthYear: The months over which the year is defined
%
% When you fix this, include the ability to select which dimension is the
% composite (just a call to permute)
%
% This function is in progress, at the moment it only a monthly range into
% yearly

% If a vector, get ts into a column vector
if isvector(ts)
    if isrow(ts)
        ts = ts';
    end
end
% ... for anything else, assume the first dimension is the progression of the
% series.

% Reshape all higher dimensional matrices into a 2D matrix
sField = size(ts);
ts = reshape( ts, sField(1), prod( sField(2:end) ));
% ... we now have a 2D matrix with a series progressing down the column

%% Choose between truncating the time series, or compositing from a certain end date

if nargin == 6
    % Get the date at which to truncate the time series
    endDate = varargin{1};
    
    % Get the index of this date
    [~, endDex, ~] = intersect(obsDatenum, endDate);
    
    % Truncate everything after the date
    if endDex < size(ts,1)
        ts( endDex+1:end,:) = [];
        obsDatenum( endDex+1:end ) = [];
        obsDatevec( endDex+1:end,: ) = [];
    end
end
    
    
    



% If there are observations before the first full year
startDex = find( obsDatevec(:,2) == monthYear(1), 1, 'first');
if startDex ~= 1
    % Truncate them from the time series
    ts( 1 : startDex-1,: ) = [];
    obsDatenum( 1:startDex-1 ) = [];
    obsDatevec( 1:startDex-1,: ) = [];
end

% Do the same for observations after the last full year
endDex = find( obsDatevec(:,2) == monthYear(2), 1, 'last');
if endDex ~= size(ts,1)
    ts( endDex+1:end,: ) = [];
    obsDatenum( endDex+1:end ) = [];
    obsDatevec( endDex+1:end,: ) = [];
end

% Check that the time series contains a full set of composite years.
% Truncate early in time if not.

% Get the number of years in the time series
[nYearsTS, nOther] = size(ts);
nYearsTS = nYearsTS / 12;

% Ensure it is modular with nYears
extraYears = mod( nYearsTS, nYears);

% Truncate the extra years from the front
if extraYears ~= 0
    ts( 1: 12*extraYears,:) = [];
    obsDatenum( 1:12*extraYears ) = [];
    obsDatevec( 1:12*extraYears,: ) = [];
end

% Get the number of composite blocks in the time series
nBlockObs = 12*nYears;
nBlocks = size(ts,1) / nBlockObs;

compTS = NaN( nBlocks, nOther );
compDatenum = NaN(nBlocks, 1);

% Get the composite means
for k = 1:nBlocks
    compTS(k,:) = mean( ts(  (k*nBlockObs-nBlockObs+1) : k*nBlockObs, :) );
    compDatenum(k) = mean( obsDatenum( (k*nBlockObs-nBlockObs+1) : k*nBlockObs) );
end
compDatevec = datevec(compDatenum);

% Reshape the data into its previous dimensional structure
compTS = reshape( compTS, [nBlocks, sField(2:end)] );

end
