%----------------------------------
% physical parameters
%----------------------------------
c = 0.5; % advection coefficient
D = 0.0; % diffusion coefficient
T = 1000; % duration of simulation
Lx = 4*2*pi; % size of the box

%----------------------------------
% numerical parameters
%----------------------------------
Nx = 128; % number of points in space
dt = 0.1; % time step for computation
dt_diag = 1.0; % time step for diagnostics
choice_derivative = 1; % centered derivatives using a stencil of 1 + 2 x choice_derivative points
% 0 = fft
% 1 = size of stencil = 1+2*1 => second order for first and second derivatives
% 2 = size of stencil = 1+2*2 => fourth order for first and second derivatives
% 3 = size of stencil = 1+2*3 =>  sixth order for first and second derivatives 
choice_scheme = 1;
% 1 = euler : explicit first  order
% 2 = rk2   : explicit second order 
% 3 = rk4   : explicit fourth order
% 4 = Lax-Wendroff : explicit ?
% 5 = ?     : implicit ?

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
