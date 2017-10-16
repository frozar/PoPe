function f = evolve_rk4(f,S,D,c,dt,dx,periodicity,choice_derivative)
	 
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
  % k3 stage

  % estimation of the system a time t0+dt/2
  f_k3 = f + dt/2 * dtf_k2;

  % derivative at time t0
  dtf_k3 = model(f_k3,S,D,c,dx,periodicity,choice_derivative);

  % ***
  % k4 stage

  % estimation of the system a time t0+dt/2
  f_k4 = f + dt * dtf_k3;

  % derivative at time t0+dt/2
  dtf_k4 = model(f_k4,S,D,c,dx,periodicity,choice_derivative);

  % ***
  % final stage

  f = f + dt * ( (1/6)*dtf_k1 + (2/6)*dtf_k2 + (2/6)*dtf_k3 + (1/6)*dtf_k4 );

