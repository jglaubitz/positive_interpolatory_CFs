%% Script to compute positive interpolatory CFs

clc, clear

%% Setting up the script 
dim = 3; % dimension: 1, 2, 3
domain = 'cube'; % domain: cube, ball, combi 
weightFun = 'C2k'; % weight function: 1, C2k, sqrt(r) 
F = 'algebraic'; % vecor space F: algebraic, trig
d = 2; % degree of exactness: 0,1,2,...
points = 'Halton'; % type of equidistributed sequence: equid, Halton 

%for d=0:10

%% Preliminary routines 
omega = generate_weightFun( weightFun, dim); % set up weight function 
[ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d ); % setup initial basis and moments 
K = length(m); % number of basis functions 

%% Compute a nonnegative LS-CF and remove zeros 
[X, w] = compute_LSCF( dim, domain, omega, basis, m, points ); % LS-CF
[X, w] = removeZeros( X, w ); % remove all zero weights and points

%% Sukzessively apply Steinitz' Austauschsatz and remove all zero weights
while K < length(w) 
    %[d, K, length(w)]
    [X, w] = Steinitz( X, w, basis ); % apply Steinitz' Austauschsatz 
    [X, w] = removeZeros( X, w ); % remove all zero weights and points
end 

%% Test 
Phi = basis(X); 
error_ex = max( abs( Phi*w-m ) );
w_min = min(w); 
[d, K, error_ex, w_min]
if w_min < 0 || error_ex >= 10^(-12)
    error('negative weight or inexact!')
end 

%% Save CF 
CF_NNI = zeros(length(w),dim+3); 
CF_NNI(:,1:dim) = X; % data points 
CF_NNI(:,dim+1) = w; % cubature weights 
CF_NNI(1,dim+2) = d; CF_NNI(1,dim+3) = K; % d and K
save( ['CFs/CF_NNI_dim=',num2str(dim),'_',domain,'_',weightFun,'_F=',F,'_d=',num2str(d),'.mat'], 'CF_NNI' ); % safe matrix

%end