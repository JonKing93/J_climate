function[Data2D, ndimSize, dimOrder] = dimNTodim2(Data, varargin)
%% Converts an n-dimensional dataset to a 2D matrix with a dimension of
% choice along the first dimension. Use this function in conjunction with
% the dim2TodimN function.
%
%
% [Data2D, ndimSize, dimOrder] = dimNTodim2(Data)
%
% [...] = dimNTodim2(Data, dim)
%
%
% ----- Inputs -----
%
% Data: An N-dimensional dataset
%
% *** Optional Inputs ***
%
% dim: (Default = 1) The dimension to set along the first dimension of the
%   2D matrix.
%
%
% ----- Output -----
%
% Data2D: The permuted data, reshaped into a 2D matrix
%
% ndimSize: The sizes of the n-dimensional, permuted dataset
%
% dimOrder: The order needed to permute the n-dimensional dataset back to
% its original ordering.

% Error check and setup
[dim, d] = setup(Data, varargin);

% Permute the data so that the dimension of interest is along the first
% dimension
dimOrder = 1:d;
dimOrder(1) = dim;
dimOrder(dim) = 1;

Data = permute(Data, dimOrder);

% Get the size of the permuted data
ndimSize = size(Data);

% Reshape the data into a 2D matrix
Data2D = reshape(Data, ndimSize(1), prod( ndimSize(2:end)));

end

% ----- Helper functions -----

function[dim,d] = setup(Data, inargs)

if length(inargs) > 1
    error('Too many input arguments');
end

% Get a value for dim
dim = 1;
if ~isempty(inargs)
    if isscalar( inargs{1})
        dim = inargs{1};
    else
        error('Unrecognized input argument');
    end
end

d = ndims(Data);
if dim > d
    error('The value for dim is greater than the number of dimensions in the dataset.');
end

end
