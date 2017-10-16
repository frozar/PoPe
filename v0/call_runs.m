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
  c_range = .25*[-2.^[-4:0] 2.^[-4:0]];%[0.0125 0.025 0.05 0.1 0.2]; % advection coefficient
  D_range = [0.0 10^-6 10^-5 10^-4 10^-3 10^-2 10^-1]; % diffusion coefficient
  T_range = 128; % duration of simulation
  Lx_range = 4*2*pi; % size of the box
  
  %----------------------------------
  % numerical parameters
  %----------------------------------
  Nx_range = 128; %[16 32 64 128 256 512 1024]; % number of points in space
  dt_range = 0.025; %[0.20 0.10 0.05 0.025 0.0125 0.00625]; % time step for computation
  dt_diag_range = 0.25; % time step for diagnostics
  choice_derivative_range = [1]; %[1 2 3 0]; % centered derivatives using a stencil of 1 + 2 x choice_derivative points
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

% choice of colors used when displaing results
colors = hsv(numel(c_range));

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
	for Nx = Nx_range
	  for dt = dt_range 
	    for dt_diag = dt_diag_range
	      for choice_derivative =  choice_derivative_range
		for choice_scheme = choice_scheme_range
  
%   %----------------------------------
%   % initial conditions
%   %----------------------------------
%   init_f = 2; % choice of the shape of f_init
%   % 1 : gaussian shape => f_init = amp_f * exp(-(x-Lx/2).^2/10);
%   % 2 : sinus          => f_init = amp_f * sin(pulse*2*pi*x/Lx+pi/2);
%   amp_f = 1;  % amplitude of f_init
%   periodicity = 1; 
%   % 0 : no periodicity
%   % 1 : periodicity
%   pulse = 3; % in case of sinus (init_f = 2)  
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


  %----------------------------------
  % "missing" parameters already determined by other parameters
  %----------------------------------
  dx = Lx/(Nx-1);
  Nt = T/dt;


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
[f_save dtf_save] = run(f_init,S,D,c,periodicity,Nt,dt,dx,dt_diag,choice_derivative,choice_scheme);
iter_run = iter_run + 1;

%----------------------------------
% 5) basic visualisation
%----------------------------------

if(1==0)

% label !
% titel !
figure(1)
subplot(1,2,1)
imagesc(f_save); colorbar;
subplot(1,2,2)
imagesc(dtf_save); colorbar
% subplot(1,3,3)
% dx2f_save = derivative2(f_save,dx,choice_derivative,periodicity);
% imagesc(dtf_save./dx2f_save); colorbar

figure(2)
fft_f_save = fft(f_save);
imagesc(abs(fft_f_save));

figure(3)
plot(sum(f_save));

end

if(1==0)
figure(4)
for i = 1:size(f_save,2)
  plot(f_init,'r--'); hold on; 
  plot(ones(size(f_init)),'g--'); hold on; 
  plot(-ones(size(f_init)),'g--'); hold on; 
  plot(S,'b--'); hold on; 
  plot(f_save(:,i),'k'); drawnow;
  pause(0.01);
  hold off;
end

end

%----------------------------------
% 6) PoPe analysis
%----------------------------------

[alpha epsilon]=PoPe_t(f_save,dtf_save,choice_model,use_source,choice_derivative_PoPe,S,periodicity,dx);

%----------------------------------
% 7) visulaisation of PoPe results
%----------------------------------

iter_color = iter_c;

figure(100);
plot(alpha(1,2:end-1),'color',colors(iter_color,:)); hold on;
plot(alpha(2,2:end-1),'--','color',colors(iter_color,:)); hold on;

tmp = mean(alpha(:,2:end-1),2)';
c_eff = tmp(1);
D_eff = tmp(2);
[c D c_eff D_eff c-c_eff D-D_eff]


figure(101);

subplot(2,2,1)
plot(c,c_eff,'*','color',colors(iter_color,:));hold on;
xlabel('c');
ylabel('c_{eff}');
plot(c_range,0*c_range,'k--');
plot(c_range,c_range,'k--');
plot(c_range,-c_range,'k--');
axis square

subplot(2,2,2)
plot(log10(abs(c)),log10(abs(c_eff)),'*','color',colors(iter_color,:));hold on;
xlabel('log10 : abs(c)');
ylabel('log10 : abs(c_{eff})');
plot(log10(c_range),log10(c_range),'k--');
axis square

subplot(2,2,3)
plot(c_eff,D_eff,'*','color',colors(iter_color,:));hold on;
xlabel('c_{eff}');
ylabel('D_{eff}');
axis square

subplot(2,2,4)
plot(log10(abs(c_eff)),log10(abs(D_eff)),'*','color',colors(iter_color,:));hold on;
xlabel('log10 : abs(c_{eff})');
ylabel('log10 : abs(D_{eff})');
axis square

figure(102);

subplot(2,2,1)
plot(D,D_eff,'*','color',colors(iter_color,:));hold on;
xlabel('D');
ylabel('D_{eff}');
plot(D_range,0*D_range,'k--');
plot(D_range,D_range,'k--');
plot(D_range,-D_range,'k--');
axis square

subplot(2,2,2)
plot(log10(abs(D)),log10(abs(D_eff)),'*','color',colors(iter_color,:));hold on;
xlabel('log10 : abs(D)');
ylabel('log10 : abs(D_{eff})');
plot(log10(D_range),log10(D_range),'k--');
axis square

subplot(2,2,3)
plot(D_eff,c_eff,'*','color',colors(iter_color,:));hold on;
xlabel('D_{eff}');
ylabel('c_{eff}');
axis square

subplot(2,2,4)
plot(log10(abs(D_eff)),log10(abs(c_eff)),'*','color',colors(iter_color,:));hold on;
xlabel('log10 : abs(D_{eff})');
ylabel('log10 : abs(c_{eff})');
axis square

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
