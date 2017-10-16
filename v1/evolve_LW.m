function f = evolve_LW(f,S,D,c,dt,dx,periodicity,choice_derivative)

  % derivative at time t0
  dtf = model(f,S,D,c,dx,periodicity,choice_derivative);

  % estimation of the system a time t0+dt
  f(2:end-1) = (f(1:end-2)+f(3:end))/2 + dt * dtf(2:end-1);
  if(periodicity==1)
    f(1)     = (f(end)+f(2))/2 + dt * dtf(1);
    f(end)   = (f(end-1)+f(1))/2 + dt * dtf(end);
  else
    f(1)     = f(1) + dt * dtf(1);
    f(end)   = f(end) + dt * dtf(end);
  end
