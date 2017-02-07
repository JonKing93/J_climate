function[grandCov] = M_SSA(Data,M)
%% Runs Multi-channel SSA
%
%
%
%
% CURRENTLY USING BROOMHEAD METHOD
%
%
%
%
%
%

[N,nchannels] = size(Data);
ntraj = N-M+1;

% Construct a grand block matrix of trajectories, here we store as a cell
% since the blocks are non-interacting.

% Preallocate
grandTraj = cell(nchannels,1);

% Fill in each channel block. Each block will be a Hankel matrix (down left
% diagonals are equal). Use this to speed construction.
for c = 1:nchannels
    grandTraj{c} = hankel( Data(1:M,c), Data(end-ntraj+1:end,c));
end


% Construct a grand covariance matrix. Here we will need multiple blocks
% for SVD, so no cell is used

% Preallocate
grandCov = NaN( nchannels*M,  nchannels*M);

% Calculate covariance blocks. Each block C_j,k is the transpose of the
% corresponding C_k,j. Use this to speed construction.
for row = 1:nchannels
    for col = row:nchannels
        % Calculate covariance
        blockCov = grandTraj{row} * grandTraj{col}' .* (1/ntraj);
        
        % Fill in the upper matrix block
        grandCov(row*M-M+1:row*M, col*M-M+1:col*M) = blockCov;
        
        if row ~= col
            % Fill in the lower matrix block with the transpose
            grandCov(col*M-M+1:col*M, row*M-M+1:row*M) = blockCov';
        end
    end
end
            
end


