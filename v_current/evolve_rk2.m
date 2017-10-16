function f = evolve_rk2(f,S,D,c,dt,dx,periodicity,choice_derivative)
	 
  % ***
  % k1 stage

  % derivative at time t0
  dtf_k1 = model(f,S,D,c,dx,periodicity,choice_derivative);

  % ***
  % k2 stage

  % estimation of the system a time t0+dt/2
  f_k2 = f + dt/2 * dtf_k1;

  % derivative at time t0+dt/2
  dtf_k2 = model(f_k2,S,D,c,dx,periodicity,choice_derivative);

  % ***
  % final stage

  f = f + dt * dtf_k2;

