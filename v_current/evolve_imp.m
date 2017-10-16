function f = evolve_imp(f,S,D,c,dt,dx,periodicity,choice_derivative,A)

  f = A\f + dt*S;
