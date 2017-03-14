function[dt] = dy2dt(dy, leapType)
%% Converts decimal year values to datetimes
%
% [dt] = dy2dt(dy, leapType)
% converts decimal year values to datetimes, including or excluding leap
% days as desired.
%
%
% -- Inputs --
%
% dy: An array of decimal year values
%
% leapType: A switch to include or exclude leap days
%   'leap': Include leap days
%   'noLeap': Exclude leap days
%

% Separate the year from the remaining time.
year = floor(dy);
day = dy-year;

% Multiply the decimal by the number of days in the year
multiplier = repmat(365, size(dy));

% If including leap years, adjust leap year multiplier to 366
if strcmpi(leapType, 'noLeap')
    % Do nothing...
elseif strcmpi(leapType, 'leap')
    % Leap years are those that are divisible by 4
    multiplier( mod(year, 4) == 0) = 366;
else
    error('Unrecognized noLeap input');
end

% Get the datetime array
dt = datetime(year,1,1) + caldays(  round(day .* multiplier)  );