function[Data] = dim2TodimN(Data2D, ndimSize, dimOrder)
%% Converts a 2D matrix to a previous N-dimesional, ordered form. Use this
% in conjunction with the dimNTodim2 function.
%
% [Data] = dim2TodimN(Data2D, ndimSize, dimOrder)
%
% ----- Inputs -----
%
% Data2D: A 2D matrix of values
%
% ndimSize: The desired size of the n-dimensional, permuted matrix. (See dimNTodim2)
%
% dimOrder: The permutation order to re-order the matrix (See dimNTodim2)
%
% ----- Output -----
%
% Data: The N-dimensional, ordered dataset

% Reshape the 2D Data
Data = reshape(Data2D, ndimSize);

% Reorder the data
Data = permute(Data, dimOrder);

end
