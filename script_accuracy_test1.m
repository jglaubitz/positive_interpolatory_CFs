%% Script to investigate accuracy 

%% Setting up the script 
clc, clear 

% free parameters
dim = 3; % dimension (1,2,3)
noise_level = 0; % 0, 10^(-6)

% fixed parameters
domain = 'cube'; % domain (cube, ball) 
weightFun = '1'; % weight function - 1, C2k, sqrt(r)
F = 'algebraic'; % vecor space F: algebraic, trig

omega = generate_weightFun( weightFun, dim); 

if dim == 1 
    f = @(x) 1./(1+x.^2); % test function 
    I = 2*atan(1); % exact integral 
elseif dim == 2 
    f = @(x,y) (1./(1+x.^2)).*(1./(1+y.^2)); % test function 
    I = (2*atan(1))^2; % exact integral 
else 
    f = @(x,y,z) (1./(1+x.^2)).*(1./(1+y.^2)).*(1./(1+z.^2)); % test function 
    I = (2*atan(1))^3; % exact integral 
end

NN_NNI = zeros(11,1); 
err_NNI = zeros(11,1);
NN_Leg = zeros(11,1); 
err_Leg = zeros(11,1);

i = 1;
for d=0:10
    
    % NNI-CF
    example = matfile(['CFs/CF_NNI_dim=',num2str(dim),'_',domain,'_',weightFun,'_F=',F,'_d=',num2str(d),'.mat']);
    C = example.CF_NNI; 
    [ M, aux] = size(C); % number of data points 
    NN_NNI(i) = M;
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = zeros(M,1); 
    for m = 1:M 
       if dim == 1  
            f_values(m) = f( X(m,1) ); 
       elseif dim == 2  
            f_values(m) = f( X(m,1), X(m,2) );
       elseif dim == 3  
            f_values(m) = f( X(m,1), X(m,2) , X(m,3) );
       else 
            error('Desired dimension not yet implemented!') 
       end
    end 
    % generate and add uniform noise 
    noise = noise_level*(2*rand(M,1)-1); 
    f_values = f_values + noise;
    CF = dot( w, f_values ); % value of the CF 
    err_NNI(i) = abs( I - CF ); % absolute error
    
    % Legendre rule 
    example = matfile(['CFs/CF_Leg_dim=',num2str(dim),'_',domain,'_d=',num2str(d),'.mat']);
    C = example.CF_Leg; 
    [ M, aux] = size(C); % number of data points 
    NN_Leg(i) = M;
    X = C(:,1:dim); % data points 
    w = C(:,dim+1); % weights 
    % function values 
    f_values = zeros(M,1); 
    for m = 1:M 
       if dim == 1  
            f_values(m) = f( X(m,1) ).*omega( X(m,1) ); 
       elseif dim == 2  
            f_values(m) = f( X(m,1), X(m,2) ).*omega( X(m,1), X(m,2) );
       elseif dim == 3  
            f_values(m) = f( X(m,1), X(m,2) , X(m,3) ).*omega( X(m,1), X(m,2) , X(m,3) );
       else 
            error('Desired dimension not yet implemented!') 
       end
    end 
    % generate and add uniform noise 
    noise = noise_level*(2*rand(M,1)-1); 
    f_values = f_values + noise;
    CF = dot( w, f_values ); % value of the CF 
    err_Leg(i) = abs( I - CF ); % absolute error
    
    i = i+1;
    
end 

    figure(1) 
    p = plot( NN_NNI,err_NNI,'r-', NN_Leg,err_Leg,'b--' ); 
    set(p, 'LineWidth',2.5)
    set(gca, 'FontSize', 18)  % Increasing ticks fontsize
    xlim([ max([NN_NNI(1);NN_Leg(1)]), min([NN_NNI(end);NN_Leg(end)]) ]) 
    %ylim([ min([err_LS2;err_l1;err_MC;err_Leg])-1, max([err_LS2;err_l1;err_MC;err_Leg])+1 ])
    xlabel('$N$','Interpreter','latex') 
    ylabel('$|C[f] - I[f]|$','Interpreter','latex')
    %set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    id = legend('interpolatory','Legendre','Interpreter','latex','Location','best');
    set(id, 'Interpreter','latex', 'FontSize',26)
    str = sprintf( ['plots/accuracy_',domain,'_dim=',num2str(dim),'_',weightFun,'.fig'] );
    savefig(str); 