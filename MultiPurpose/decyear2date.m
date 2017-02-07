function[timenum, timevec] = decyear2date(yeardec, varargin)
%% This function converts decimal time to a datenum. The function may be
% specified to include or skip a leap year.

if ~isvector(yeardec)
    error('Yeardec must be a vector');
end
if ~iscolumn(yeardec)
    yeardec = yeardec';
end
    

% Default is to include leap year
if nargin == 1
    includeLeap = true;
    
% Otherwise check user input for no leap year.
elseif nargin == 2
    if strcmpi(varargin{1}, 'noleap')
        includeLeap = false;
    elseif strcmpi( varargin{1}, 'leap')
        includeLeap = true;
    else
        error('Unrecognized input');
    end
else
    error('Too many inputs');
end

% Separate the year and decimal
year = fix(yeardec);
day = yeardec - year;

% Get the number of days as a multiplier
nTime = length(yeardec);
multiplier = repmat(365, nTime, 1);

% If including leap years, make the leap year multiplier 366
if includeLeap
    multiplier( mod(year,4) == 0) = 366;
end
    
% Convert to a datenum
timenum = datenum(year,1,1) + (day .* multiplier);
timevec = datevec(timenum);
end

