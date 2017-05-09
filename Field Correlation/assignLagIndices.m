function[i1, i2] = assignLagIndices(dates1, dates2)
%% Assigns time indices for lagFieldcorr to two sets of overlapping, 
% equally spaced datetime vectors. The first overlapping date will be given
% an index of 0.
%
% [i1, i2] = assignLagIndices(dates1, dates2)
%
%
% ----- Inputs -----
%
% dates1: A datetime vector. Must be equally spaced.
%
% dates2: A datetime vector. Must be equally spaced at the same interval as
%       dates1 and overlap dates1.
%
%
% ----- Outputs -----
%
% i1: The lag indices for the first set of dates.
%
% i2: The lag indices for the second set of dates.
%
%
% ----- Written By -----
%
% Jonathan King, 2017, University of Arizona, jonking93@email.arizona.edu

% Initial error checking
errCheck(dates1, dates2);

% Check if the beginning of the first series overlaps the second
if any( dates2 == dates1(1) )
    
    % First overlap is zero, and occurs at beginning of dates1
    i1 = ( 0:length(dates1)-1 )';
    
    % Get the index of overlap for dates2
    dex = find(dates2 == dates1(1));
    
    % Get the indices for dates2
    i2 = (  -(dex-1):length(dates2)-dex  )';
    
    
else
    
    % First overlap is at beginning of dates2
    i2 = ( 0:length(dates2)-1 )';
    
    % Get the index of overlap for dates2
    dex = find(dates1 == dates2(1));
    
    % Get the indices for dates1
    i1 = (  -(dex-1):length(dates1)-dex  )';
end

end

% ----- Helper Method -----

function[] = errCheck(dates1, dates2)

if ~isvector(dates1) || ~isvector(dates2) || ~isdatetime(dates1) || ~isdatetime(dates2)
    error('dates1 and dates2 must be datetime vectors');

elseif hasNaN(dates1) || hasNaN(dates2)
    error('dates1 and dates2 may not contain NaN');

else
    
    spacing1 = dates1(2:end) - dates1(1:end-1);
    spacing2 = dates2(2:end) - dates2(1:end-1);
    
    if any(spacing1(2:end)~=spacing1(1:end-1))
        error('dates1 does not have a constant time interval')
    elseif any(spacing2(2:end)~=spacing2(1:end-1))
        error('dates2 does not have a constant time interval')
    elseif spacing1(1) ~= spacing2(1)
        error('dates1 and dates2 must have the same time interval');
    end
end
end

