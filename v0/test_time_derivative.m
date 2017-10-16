% this script tests if the time derivative is properly computed and store in the run.m script

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
  init_f = 2; % choice of the shape of f_init
  % 1 : gaussian shape => f_init = amp_f * exp(-(x-Lx/2).^2/10);
  % 2 : sinus          => f_init = amp_f * sin(pulse*2*pi*x/Lx+pi/2);
  amp_f = 1;  % amplitude of f_init
  periodicity = 1; 
  % 0 : no periodicity
  % 1 : periodicity
  pulse = 3; % in case of sinus (init_f = 2)  
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
[x f_init S] = init(Nx,Lx,init_f,amp_f,pulse,init_S,amp_S);

%----------------------------------
% 4) run simulation
%----------------------------------
step = dt_diag/dt;
[f_save_1 dtf_save_1 time_1] = run(f_init,S,D,c,periodicity,Nt,dt,dx,dt_diag,choice_derivative,choice_scheme);
dt_diag = dt;
[f_save_2 dtf_save_2 time_2] = run(f_init,S,D,c,periodicity,Nt,dt,dx,dt_diag,choice_derivative,choice_scheme);


% compute dtf offline

dtf_save_2(:,1+2:end-2) = (1*f_save_2(:,1+2-2:end-2-2) - 8*f_save_2(:,1+2-1:end-2-1) +8*f_save_2(:,1+2+1:end-2+1) -1*f_save_2(:,1+2+2:end-2+2))/(12*dt);

error_dtf_save = dtf_save_2(:,1:step:end-step) - dtf_save_1(:,1:end-1);

figure(1)
subplot(3,1,1);
imagesc(dtf_save_1(:,1:end-1)); colorbar;
subplot(3,1,2);
imagesc(dtf_save_2(:,1:step:end-step)); colorbar;
subplot(3,1,3);
imagesc(error_dtf_save); colorbar
