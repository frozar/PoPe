% clean everything in workspace
close all; clear all; clc;


%----------------------------------
% 1) select a special dataset for the initialisation of the system
%----------------------------------
INIT = 0;
% 0 = custom : everything has to be specified in the switch over INIT
% 1 = FTCS_CRASH : inputs are read from a the file data_FTCS_CRASH
%     euler time integration + centered finite differences for transport equation
%     -> crash of the simulation
% 2 =

switch INIT

case 0 % custom

  %----------------------------------
  % physical parameters
  %----------------------------------
  c_range = 0;%
  D_range = 10^-2; % diffusion coefficient
  T_range = 16; % duration of simulation
  Lx_range = 4*2*pi; % size of the box
  
  %----------------------------------
  % numerical parameters
  %----------------------------------
  Nx_range = [256]; % number of points in space
  dt_range = 0.01*2.^[0:-1:-3]; % time step for computation
  dt_diag_range = 0.1; % time step for diagnostics
  choice_derivative_range = [1 2 3]; %[1 2 3 0]; % centered derivatives using a stencil of 1 + 2 x choice_derivative points
  % 0 = fft
  % 1 = size of stencil = 1+2*1 => second order for first and second derivatives
  % 2 = size of stencil = 1+2*2 => fourth order for first and second derivatives
  % 3 = size of stencil = 1+2*3 =>  sixth order for first and second derivatives 
  choice_scheme_range = [1];
  % 1 = euler : explicit first  order
  % 2 = rk2   : explicit second order 
  % 3 = rk4   : explicit fourth order
  % 4 = Lax-Wendroff : explicit
  % 5 = ?     : implicit order 1
    
  %----------------------------------
  % initial conditions
  %----------------------------------
  %%% initial state of the system
  init_f = 2; % choice of the shape of f_init
  % 1 : gaussian shape => f_init = amp_f * exp(-(x-Lx/2).^2/10);
  % 2 : sinus          => f_init = amp_f * sin(pulse*2*pi*x/Lx+pi/2);
  % 3 : randn, normal random variable
  amp_f = 1;  % amplitude of f_init
  periodicity = 1; 
  % 0 : no periodicity
  % 1 : periodicity
  pulse = 3; % in case of sinus (init_f = 2)  
  %%% source 
  init_S = 3; %
  % 1 : gaussian shape centered in Lx/4   => S = amp_S * exp(-(x-Lx/4).^2/10);
  % 2 : gaussian shape centered in 3*Lx/4 => S = amp_S * exp(-(x-3*Lx/4).^2/10);
  % 3 : sum of shape #1 and shape #2
  amp_S = 0;

  %----------------------------------
  % paramaters for the PoPe analysis
  %----------------------------------
  
  choice_model_range = [3];
  % 1 : only dissipation operator
  % 2 : only advection operator
  % 3 : dissipation and advection operators
  use_source_range = 0;
  % 0 : does not use any source
  % 1 : use the source
  
  choice_derivative_PoPe_range = 0; % centered derivatives using a stencil of 1 + 2 x choice_derivative points
  % 0 = fft
  % 1 = size of stencil = 1+2*1 => second order for first and second derivatives
  % 2 = size of stencil = 1+2*2 => fourth order for first and second derivatives
  % 3 = size of stencil = 1+2*3 =>  sixth order for first and second derivatives 
  
case 1 % use file FTCS_CRASH :  
     data_FTCS_CRASH

case 2 % ?

end

% structure containing results of each run : output of simulations and outputs of PoPe treatment
data = struct('f_save',[],'dtf_save',[],'time',[],'x',[],'nb_alpha',[],'alpha_th',[],'alpha',[],'epsilon',[],'dt',[],'dx',[]);

%----------------------------------
% nested loop to explore parameters
%----------------------------------
iter_run = 0; % iterator over the runs

%----------------------------------
% physical parameters
%----------------------------------
iter_c = 0; % iterator over the c-values
for c = c_range
  iter_c = iter_c + 1;
  for D = D_range
    for T = T_range
      for Lx = Lx_range
	  
	%----------------------------------
	% numerical parameters
	%----------------------------------
	iter_Nx = 0; % iterator over the Nx-values
	for Nx = Nx_range
	  iter_Nx = iter_Nx +1;
	  iter_dt = 0; % iterator over the dt-values
	  for dt = dt_range 
	    iter_dt = iter_dt +1;
	    %dt = dt/Nx;
	    for dt_diag = dt_diag_range
	      for choice_derivative =  choice_derivative_range
		for choice_scheme = choice_scheme_range
  
%   %----------------------------------
%   % initial conditions
%   %----------------------------------
%   %%% initial state of the system
%   init_f = 2; % choice of the shape of f_init
%   % 1 : gaussian shape => f_init = amp_f * exp(-(x-Lx/2).^2/10);
%   % 2 : sinus          => f_init = amp_f * sin(pulse*2*pi*x/Lx+pi/2);
%   amp_f = 1;  % amplitude of f_init
%   periodicity = 1; 
%   % 0 : no periodicity
%   % 1 : periodicity
%   pulse = 3; % in case of sinus (init_f = 2)  
%   %%% source 
%   init_S = 3; %
%   % 1 : gaussian shape centered in Lx/4   => S = amp_S * exp(-(x-Lx/4).^2/10);
%   % 2 : gaussian shape centered in 3*Lx/4 => S = amp_S * exp(-(x-3*Lx/4).^2/10);
%   % 3 : sum of shape #1 and shape #2
%   amp_S = 0;

		  %----------------------------------
		  % paramaters for the PoPe analysis
		  %----------------------------------
		  for choice_model = choice_model_range
		    for use_source = use_source_range  
		      for choice_derivative_PoPe = choice_derivative_PoPe_range

% iterator over the different runs
iter_run = iter_run + 1;

% save the theoretical value of weigths : alpha_th
switch choice_model
  case 1  % only dissipation operator
       alpha_th = D;
  case 2  % only advection operator
       alpha_th = c;
  case 3  % dissipation and advection operators
       alpha_th = [c D];
end
if(use_source==1)
  alpha_th = [alpha_th 1];
end
data(iter_run).alpha_th=alpha_th;
data(iter_run).nb_alpha = numel(alpha_th);

  %----------------------------------
  % "missing" parameters already determined by other parameters
  %----------------------------------
  dx = Lx/(Nx-1);
  Nt = T/dt;

  % save important parameters
  data(iter_run).dt = dt;
  data(iter_run).dx = dx;

%----------------------------------
% 2) verification of the theoretical stability
%----------------------------------
% CFL !

%----------------------------------
% 3) initialisation of the system
%----------------------------------
[x f_init S] = init_system(Nx,Lx,init_f,amp_f,pulse,init_S,amp_S);

%----------------------------------
% 4) run simulation
%----------------------------------
[f_save dtf_save time] = run(f_init,S,D,c,periodicity,Nt,dt,dx,dt_diag,choice_derivative,choice_scheme);
data(iter_run).f_save = f_save;
data(iter_run).dtf_save = dtf_save;
data(iter_run).time = time;
data(iter_run).x = x;

%----------------------------------
% 5) PoPe analysis
%----------------------------------
[alpha epsilon]=PoPe_study_over_time(f_save,dtf_save,choice_model,use_source,choice_derivative_PoPe,S,periodicity,dx);
data(iter_run).alpha=alpha;
data(iter_run).epsilon=epsilon;

		      end
		    end
		  end
		end
	      end
	    end
	  end
	end
      end
    end
  end
end

% choice of colors used when displaing results
colors_1 = hsv(numel(data)); % first ensemble of colors
colors_2 = [1 0 0; % second ensemble of colors
	    0 1 0;
	    0 0 1;
	    0 0 0;];

%----------------------------------
% 6) basic visualisation
%----------------------------------
% basic_visualisation;

%----------------------------------
% 7) visualisation of PoPe results (over all runs)
%----------------------------------
PoPe_visualisation;
