function[sM] = sampleMonthly(date, dayType, leapType)
%% Samples a set of times on monthly intervals.
%
% [sM] = sampleMonthly(time, timeType, dayType, includeLeap)
% samples a time vector, which has a type specified by timeType, on monthly
% intervals at either the first, last, or middle day of each month, as
% specified by dayType. Adjusts dates for leap vs no leap day years.
%
% -- In --
%
% date: A datetime array
%
% dayType: Specifies the beginning, middle, or end of a month
%   'first': Sample on the first day of each month
%   'mid': Sample on the middle day of the month
%   'last': Sample on the last day of the month
%
% leapType: Specifies whether to consider leap days
%   'leap': Includes leap days
%   'noLeap': Excludes leap days
%
% -- Out --
%
% sM: The datetimes for monthly sampling points.

% Get the earliest and latest year
startYear = year( min(date));
endYear = year( max(date));

% Preallocate the year, month, and day values for the final datetime
dYear = NaN( (endYear-startYear+1)*12,1); % preallocate
dMonth = dYear;
dDay = dYear;

% Fill the values
index = 1;
for k = startYear:endYear
    dYear(index:index+11) = k;
    dMonth(index:index+11) = (1:12)';
    
    % Choose which day to assign the date
    if strcmpi(dayType,'first')
        dDay(index:index+11) = 1;
    elseif strcmpi(dayType, 'mid')
        dDay(index:index+11) = [16 15 16 16 16 16 16 16 16 16 16 16]';
    elseif strcmpi(dayType, 'last')
        dDay(index:index+11) = [31 28 31 30 31 30 31 31 30 31 30 31]';
        
        % Adjust for leap days if desired
        if strcmpi(leapType, 'leap')
            % Do nothing...
        elseif strcmpi(leapType, 'noLeap')
            dDay(index+1) = 29;
        else
            error('Unrecognized leapType');
        end
    else
        error('Unrecognized dayType');
    end
    
    % Update the indices to the next year
    index = index+12;
end

% Make the datetime array
sM = datetime(dYear, dMonth, dDay);

% Remove values that exceed the original set of dates
sM( sM < min(date) ) = [];
sM( sM > max(date) ) = [];





