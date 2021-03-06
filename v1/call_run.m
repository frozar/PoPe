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
  c = 0.1; % advection coefficient
  D = 0.0; % diffusion coefficient
  T = 100; % duration of simulation
  Lx = 4*2*pi; % size of the box
  
  %----------------------------------
  % numerical parameters
  %----------------------------------
  Nx = 128; % number of points in space
  dt = 0.1; % time step for computation
  dt_diag = 1.0; % time step for diagnostics
  choice_derivative = 0; % centered derivatives using a stencil of 1 + 2 x choice_derivative points
  % 0 = fft
  % 1 = size of stencil = 1+2*1 => second order for first and second derivatives
  % 2 = size of stencil = 1+2*2 => fourth order for first and second derivatives
  % 3 = size of stencil = 1+2*3 =>  sixth order for first and second derivatives 
  choice_scheme = 5;
  % 1 = euler : explicit first  order
  % 2 = rk2   : explicit second order 
  % 3 = rk4   : explicit fourth order
  % 4 = Lax-Wendroff : explicit
  % 5 = ?     : implicit order 1
  
  %----------------------------------
  % "missing" parameters already determined by other parameters
  %----------------------------------
  dx = Lx/(Nx-1);
  Nt = T/dt;
  
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
  
  choice_model = 3;
  % 1 : only dissipation operator
  % 2 : only advection operator
  % 3 : dissipation and advection operators
  use_source = 0;
  % 0 : does not use any source
  % 1 : use the source
  
  choice_derivative_PoPe = 0; % centered derivatives using a stencil of 1 + 2 x choice_derivative points
  % 0 = fft
  % 1 = size of stencil = 1+2*1 => second order for first and second derivatives
  % 2 = size of stencil = 1+2*2 => fourth order for first and second derivatives
  % 3 = size of stencil = 1+2*3 =>  sixth order for first and second derivatives 
  
case 1 % use file FTCS_CRASH :  
     data_FTCS_CRASH

case 2 % ?

end

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

%----------------------------------
% 5) basic visualisation
%----------------------------------
basic_visualisation;

%----------------------------------
% 6) PoPe analysis
%----------------------------------
% set the theoretical value of weigths : w_th
switch choice_model
  case 1  % only dissipation operator
       w_th = D;
  case 2  % only advection operator
       w_th = c;
  case 3  % dissipation and advection operators
       w_th = [c D];
end
if(use_source==1)
  w_th = [w_th 1];
end
%
[alpha epsilon]=PoPe_study_over_time(f_save,dtf_save,choice_model,use_source,choice_derivative_PoPe,S,periodicity,dx);
% study of w
PoPe_c(iter_run,:) = mean(alpha(:,2:end-1),2);
PoPe_Delta_C(iter_run,:) = abs(w_th - PoPe_c(iter_run,:));
PoPe_delta_c(iter_run,:) = mean(abs(PoPe_c(iter_run,:)'*ones(1,size(alpha(:,2:end-1),2))-alpha(:,2:end-1)),2);
dx_run(iter_run) = dx;
dt_run(iter_run) = dt;

%----------------------------------
% 7) visulaisation of PoPe results
%----------------------------------
iter_run = 1;
PoPe_visualisation;
