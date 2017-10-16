function f = evolve_euler(f,S,D,c,dt,dx,periodicity,choice_derivative)

  % derivative at time t0
  dtf = model(f,S,D,c,dx,periodicity,choice_derivative);

  % estimation of the system a time t0+dt
  f = f + dt * dtf;

