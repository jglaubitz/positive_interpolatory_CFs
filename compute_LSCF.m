%% compute_LSCF
% Compute a nonnegative LS-CF
% 
% INPUT: 
%  Sample : Sample for the data points 
%  A_init : Matrix containing the values of the initial basis 
%  m_init : Vector containing the moments of the initial basis
%  d :      Maximal degree 
%
% OUTPUT: 
%  X : Matrix which contains the data points 
%  w : Vector of cubature weights

function [ X, w] = compute_LSCF( dim, domain, omega, basis, m_init, points )

    K = length(m_init); % number of basis functions 

    %% routine to determine a nonnegative LS-CF 
    M = K; w_min = -1;
    while w_min < 0 
        %% Determine data points and initial Vandermonde matrix
        Sample = generate_points( points, domain, dim, omega, M ); % data points and discrete weights
        R = diag(Sample.r); % diagonal weight matrix 
        Phi = basis(Sample.coord); % initial Vandermonde matrix 
        %% Compute the LS weights 
        if K <= Sample.N && rank(Phi) == K 
            % Compute the V-matrix and moments of the DOBs 
            [ Phi, m] = GramSchmidt( Sample, Phi, m_init ); % orthonormalization 
            [ Phi, m] = GramSchmidt( Sample, Phi, m ); % re-orthonormalization 
            % ON_error = max(abs( Phi*R*(Phi') - eye(K,K) ), [], 'all') 
            % Compute the LS weights 
            w = R*(Phi')*m; % LS weights 
            w_min = min(w); % their smallest value 
            %ex_error = max( abs( Phi*w - m ) ) 
        end
        M = 2*M; % double the data points
    end

    X = Sample.coord; % data points 
    
end