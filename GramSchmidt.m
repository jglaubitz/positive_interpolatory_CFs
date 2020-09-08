%% GramSchmidt
% Generates a DOB \pi_k (evaluates them at the data points) and their moments. 
% We use the modified Gram-Schmidt process! 
% 
% INPUT: 
%  Sample : Sample for the data points 
%  A_init : Matrix containing the values of the initial basis 
%  m_init : Vector containing the moments of the initial basis
%  d :      Maximal degree 
%
% OUTPUT: 
%  A : Matrix which contains the values of \pi_k at the data points 
%  m : Vector of moments 

function [ A, m] = GramSchmidt( Sample, A_init, m_init )

    [K, N] = size(A_init); % number of basis functions and data points 
    A = A_init; % initiate the Vandermonde matrix 
    m = m_init; % initiate the vector of moments 

    % Modified Gram-Schmidt procedure 
    for k=1:K 
        for l=1:k-1 
            inner_prod = dot( A(k,:).*A(l,:), Sample.r ); % discrete inner prod
            A(k,:) = A(k,:) - inner_prod*A(l,:); % orthogonalization 
            m(k) = m(k) - inner_prod*m(l); % also change the moments 
        end
        discr_norm = sqrt( dot( A(k,:).^2 , Sample.r ) ); % discrete norm 
        A(k,:) = A(k,:)/discr_norm; % normalization 
        m(k) = m(k)/discr_norm; % also change the moments
    end
    
end