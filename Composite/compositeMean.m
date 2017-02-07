function[] = compositeMean(data, window, varargin)
%% This function gets the composite (discrete average) value of an n-dimensional 
% dataset along a dimension of choice.
%
%
%
% ----- Inputs -----
%
% data: An n-dimensional dataset
%
% window: The number of indices over which to calculate each composite
% mean.
%
% *** Optional Inputs ***
% 
% dim: (Default = 1) The dimension along which to take the mean
%
% ----- Outputs -----

%% Initial error checking and sizing
[ldata, sdata, ddata, dim] = setup(data, window, varargin);


% Move the dimension of interest to the first dimension
ddata(1) = dim;
ddata(dim) = 1;

% Permute the data matrix accordingly
data = permute( data, ddata);

% Reshape into a 2D matrix


% Preallocate the composite mean
nComps = ldata / window;
sdata(dim) = nComps;
compMean = NaN(sdata);

% Take the composite means
for k = 1:nComps
    compMean(:, 


