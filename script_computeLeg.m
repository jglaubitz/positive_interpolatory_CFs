%% Script to compute (transformed) product Legendre rules 

%% Setting up the script 
clc, clear 

dim = 3; % dimension (1,2,3)
domain = 'cube'; % domain (cube, ball) 
weightFun = 'C2k'; % weight function: 1, C2k, sqrt(r) 
F = 'algebraic'; % vecor space F: algebraic, trig


for d=0:10
    
    % load NNI-CF to use approximately the same number of data points 
    example = matfile(['CFs/CF_NNI_dim=',num2str(dim),'_',domain,'_',weightFun,'_F=',F,'_d=',num2str(d),'.mat']);
    C = example.CF_NNI; 
    [ M, aux] = size(C);
    n = ceil(M^(1/dim));
    
    % compute the Legendre rule 
    [X, w_Leg, d_Leg, K_Leg ] = compute_LegendreRule( dim, domain, n );
    
    % save points, weights, and d and K in a matrix
    CF_Leg = zeros(n^dim,dim+3); % initiate 
    CF_Leg(:,1:dim) = X; % store data points
    CF_Leg(:,dim+1) = w_Leg; % store cubature weights 
    CF_Leg(1,dim+2) = d_Leg; % store d
    CF_Leg(2,dim+3) = K_Leg; % store K      
    save( ['CFs/CF_Leg_dim=',num2str(dim),'_',domain,'_d=',num2str(d),'.mat'], 'CF_Leg' ); % safe matrix   
    
end