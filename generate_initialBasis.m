%% generate_initialBasis
% Generates the initial basis and the corresponding moments 
% 
% INPUT: 
%  dim :        dimension 
%  domain :     domain 
%  weightFun :  weight function 
%  F :          vector space F 
%  d :          degree of exactness 
%
% OUTPUT: 
%  basis : vector-valued function with basis elements 
%  m :     corresponding moments 

function [ basis, m ] = generate_initialBasis( dim, domain, weightFun, F, d )

    if strcmp( F, 'algebraic')  
        [ basis, m ] = generate_initialBasis_algPoly( dim, domain, weightFun, d ); 
    elseif strcmp( F, 'trig')
        [ basis, m ] = generate_initialBasis_trigPoly( dim, domain, weightFun, d );
    else
        error('Desired vector space F not yet implemented!')
    end
    
end