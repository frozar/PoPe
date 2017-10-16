function [f_save dtf_save time] = run(f_init,S,D,c,periodicity,Nt,dt,dx,dt_diag,choice_derivative,choice_scheme)

% main variable set to initial condition
f = f_init;
% variable for the computation of the time derivative used in PoPe
dtf = zeros(size(f_init));
% allocation of memory to save a series of f
f_save = zeros(numel(f_init),Nt*dt/dt_diag+1);
% allocation of memory to save a series of d/dt(f)
dtf_save = zeros(numel(f_init),Nt*dt/dt_diag+1);

iter_save = 1; % iterator
T = 0; % current time
time = zeros(1,Nt*dt/dt_diag+1); % time serie of diagnostics
time(:,iter_save)=T; % save the initial "current time"
f_save(:,iter_save)=f_init; % save the initial condition

% precomputation of matrix used in implicite method
if(choice_scheme==5)

  % periodicity

  A = init_evolve_imp(numel(f_init),dx,dt,c,D,periodicity);

%  figure(1000000)
%  subplot(1,2,1)
%  spy(A);axis square;
%  subplot(1,2,2)
%  imagesc(log10(abs(A)));axis square;

%  drawnow;
%  pause

end


tic
for i = 1:Nt

    %----------------------------------
    % choice of time integration scheme
    %----------------------------------
    switch choice_scheme
      case 1 % euler
	f = evolve_euler(f,S,D,c,dt,dx,periodicity,choice_derivative);
      case 2 % rk2
	f = evolve_rk2(f,S,D,c,dt,dx,periodicity,choice_derivative);
      case 3 % rk4
	f = evolve_rk4(f,S,D,c,dt,dx,periodicity,choice_derivative);
      case 4 % Lax Wendroff
	f = evolve_LW(f,S,D,c,dt,dx,periodicity,choice_derivative);
      case 5 % implicite magic
	f = evolve_imp(f,S,D,c,dt,dx,periodicity,choice_derivative,A);
    end

    % current time
    T = i*dt;

    %----------------------------------
    % save f
    %----------------------------------
    if(mod(T,dt_diag)==0)
      iter_save = iter_save +1;
      f_save(:,iter_save)=f;      
      time(:,iter_save)=T;
%      if(mod(iter_save,(Nt*dt/dt_diag)/10)==0)
%	display(['iteration = ' num2str(i) '/' num2str(Nt) ', time = ' num2str(toc)])
%      end
    end

    %----------------------------------
    % save dt f for PoPe post treatment
    %----------------------------------
    tol = dt/100; % in case of rounding error
    if(mod(T+2*dt,dt_diag)<=tol)
      % t = dt_diag - 2 dt
      dtf =      +1 * f;
    elseif(mod(T+1*dt,dt_diag)<=tol)
      % t = dt_diag - 1 dt
      dtf = dtf - 8 * f;
    %elseif(mod(T+0*dt,dt_diag)<=tol)
      % t = dt_diag
      % dtf = dtf + 0 * f;
    elseif(mod(T-1*dt,dt_diag)<=tol)
    % t = dt_diag + 1 dt
      dtf = dtf + 8 * f;
    elseif(mod(T-2*dt,dt_diag)<=tol)
      % t = dt_diag + 2 dt
      dtf = dtf - 1 * f;
      %
      dtf_save(:,iter_save)=dtf/(12*dt);
    %
    end
    
end

% time derivative of the first iteration and last iteration is not computed :
% it's slitly more difficult and not nessecary
dtf_save(:,1)=0; % initial condition / first iteration
dtf_save(:,end)=0; % last iteration
