function dtf = model(f,S,D,c,dx,periodicity,choice_derivative)
  % compute the equation of diffusion + advection with two sources
  % d/dt(f) = c d/dx(f) + D d2/dx2(f) + S
	 
%  if (D==0)
%    % just an advection equation
%    dtf = c * derivative1(f,dx,choice_derivative,periodicity) + S; 
%  elseif (c==0)
%    % just a diffusion equation	 
%    dtf = D * derivative2(f,dx,choice_derivative,periodicity) + S; 
%  else
    % a diffusion + advection equation
    dtf = c * derivative1(f,dx,choice_derivative,periodicity) + D * derivative2(f,dx,choice_derivative,periodicity) + S; 
%  end
  
    if (periodicity==0)
      dtf(1)=0;
      dtf(end)=0;
    end
