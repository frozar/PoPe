function [alpha epsilon]=PoPe_study_over_time(f,dtf,choice_model,use_source,choice_derivative_PoPe,S,periodicity,dx)

[Nx Nt]=size(f);

% model containing
switch choice_model
  case {1,2}
    % only dissipation operator
    % or only advection operator
    % thus only one weight
    nb_w = 1;
  case 3
    % dissipation and advection operators
    % thus two weigths
    nb_w = 2;
end
% add the number of sources taking into account
nb_w = nb_w + use_source;

delta_x = 5; % points not considered in the analysis at boundaries in space

% allocate memory
alpha = zeros(nb_w,Nt); % effective weights of operators
epsilon = zeros(Nx,Nt); % discrepency between effective equations and theoretical equation
Operators = zeros(nb_w,Nx);

for i = 2:Nt-1

  iter_operators = 0;

  % post computation of operators
  switch choice_model
    case 1
      % only dissipation operator
      iter_operators = iter_operators + 1; 
      Operators(iter_operators,:) = derivative2(f(:,i),dx,choice_derivative_PoPe,periodicity);
    case 2
      % only advection operator
      iter_operators = iter_operators + 1; 
      Operators(iter_operators,:) = derivative1(f(:,i),dx,choice_derivative_PoPe,periodicity);
    case 3
      % dissipation and advection operators
      % thus two weigths
      iter_operators = iter_operators + 1; 
      Operators(iter_operators,:) = derivative1(f(:,i),dx,choice_derivative_PoPe,periodicity);
      iter_operators = iter_operators + 1; 
      Operators(iter_operators,:) = derivative2(f(:,i),dx,choice_derivative_PoPe,periodicity);
  end
  
  % taking into account the source
  if(use_source==1)
    iter_operators = iter_operators + 1;   
    Operators(iter_operators,:) = S;
  end


  if(1==1)

    % matrix in least mean square approach  
    A = Operators(:,1+delta_x:end-delta_x) * Operators(:,1+delta_x:end-delta_x)';
    RHS = Operators(:,1+delta_x:end-delta_x) * dtf( 1+delta_x : end-delta_x , i );
    % projection in least mean square approach
    alpha(:,i) = A\RHS;
    epsilon(:,i) = dtf(:,i)' - alpha(:,i)'*Operators;

  else

    % matrix
    A = Operators(:,1+delta_x:end-delta_x)';
    RHS = dtf( 1+delta_x : end-delta_x , i );
    % projection
    alpha(:,i) = A\RHS;
    epsilon(:,i) = dtf(:,i)' - alpha(:,i)'*Operators;
    
%    size(A)
%    size(RHS)
%    size(A\RHS)

  end
   
end
