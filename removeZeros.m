%% removeZeros
% Compute a nonnegative LS-CF
% 
% INPUT: 
%  X : Matrix which contains the data points 
%  w : Vector of cubature weights
%
% OUTPUT: 
%  Y : Matrix which contains the data points 
%  v : Vector of cubature weights (no zeros anymore)

function [ X, w] = removeZeros( X, w )

    [N,dim] = size(X); 
    n = 1; 
    while n <= N 
        if abs( w(n) ) <= 10^(-15) % zero weight 
            X(n,:) = []; % remove the data points
            w(n) = []; % remove the weight 
            N = N-1; 
        else 
            n = n+1; 
        end
    end 
    
end