function [x f_init S] = init(Nx,Lx,init_f,amp_f,pulse,init_S,amp_S)

%----------------------------------
% initial condition of the main variable : f
%----------------------------------
switch init_f

  case 1 % a gaussian shape
    x = [0:Lx/(Nx-1):Lx]';
    f_init = amp_f * exp(-(x-Lx/2).^2/10);
    
  case 2 % a sinus
    nx = Nx+1;
    x_tmp = [0:Lx/(nx-1):Lx]';
    x = x_tmp(1:end-1);
    f_init = amp_f * sin(pulse*2*pi*x/Lx+pi/2);

  case 3 % a sinus
    x = [0:Lx/(Nx-1):Lx]';
    f_init = rand(size(x));

end

%----------------------------------
% define the source : S
%----------------------------------
switch init_S

  case 1
    S = amp_S * exp(-(x-Lx/4).^2/10);

  case 2
    S = amp_S * exp(-(x-3*Lx/4).^2/10);

  case 3
    S = amp_S * (exp(-(x-Lx/4).^2/10) - exp(-(x-3*Lx/4).^2/10));

end




