function[timenum, timevec] = monthlyDate(startYear, endYear, dayType, edgeMonths, includeLeap)
%% Gets the date for either the first day, middle, or last day of the
% months from a starting year to an ending year. 
%
% [timenum, timevec] = monthlyDate(startYear, endYear, monthType, edgeMonths, includeLeap)%

if strcmpi(includeLeap, 'leap')
    includeLeap = true;
else 
    includeLeap = false;
end

% Initialize a vector of the sampled years
nYears = endYear - startYear + 1;
years =  NaN( 12*nYears, 1);
for k = 1:nYears
    years( k*12-11:k*12 ) = startYear + k - 1;
end

% Get a vector of the sampled months
months = repmat( [1:12]', nYears, 1); %#ok<NBRAK>

% Get a vector of the sampled day
switch dayType
    
    % The first day of each month
    case 'first'
        days = ones( 12*nYears, 1);
        
    % The middle day of each month (rounds up for even numbered months)
    case 'mid'
        days = repmat( [16 15 16 16 16 16 16 16 16 16 16 16]', nYears, 1);
        
    % The last day of each month    
    case 'last'
        days = repmat( [31 28 31 30 31 30 31 31 30 31 30 31]', nYears, 1);
        
        % Switch February to 29 if including leap years
        if includeLeap
            days(  mod(years, 4)==0  &  months==2  ) = 29;
        end
        
    otherwise
        error('Unrecognized day type');
end

% Convert to datenums
timenum = datenum(years, months, days);

% Eliminiate the months outside of the edge months
if edgeMonths(1) > 1
    timenum( 1 : edgeMonths(1)-1 ) = [];
end
if edgeMonths(2) < 12
    timenum(   end- (12 - edgeMonths(2)) +1 : end) = [];
end

% Also get a datevec
timevec = datevec(timenum);
        
end
        
        
        