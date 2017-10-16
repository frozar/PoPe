function A = init_evolve_imp(Nx,dx,dt,c,D,periodicity)

  if(periodicity==0)
    
    D1_nonperiodic = (diag(ones(1,Nx-1),+1) - diag(ones(1,Nx-1),-1))/(dx*2);
    D1_nonperiodic(1,:) = 0;
    D1_nonperiodic(end,:) = 0;
    
    D2_nonperiodic = (diag(ones(1,Nx-1),+1) - 2 * diag(ones(1,Nx)) + diag(ones(1,Nx-1),-1))/(dx^2);
    D2_nonperiodic(1,:) = 0;
    D2_nonperiodic(end,:) = 0;
    
    A_nonperiodic = diag(ones(1,Nx)) - dt * ( c*D1_nonperiodic + D*D2_nonperiodic);
    A = A_nonperiodic;

  else
      
    D1_periodic    = (diag(ones(1,Nx-1),+1) - diag(ones(1,Nx-1),-1))/(dx*2);
    D1_periodic(1,end) = -1/(dx*2);
    D1_periodic(end,1) = 1/(dx*2);
    
    D2_periodic    = (diag(ones(1,Nx-1),+1) - 2 * diag(ones(1,Nx)) + diag(ones(1,Nx-1),-1))/(dx^2);
    D2_periodic(1,end) = 1/(dx^2);
    D2_periodic(end,1) = 1/(dx^2);
    
    A_periodic    = diag(ones(1,Nx)) - dt * ( c*D1_periodic    + D*D2_periodic);
    A = A_periodic;    

  end
