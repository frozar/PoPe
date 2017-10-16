close all; clear all; clc;
format long

%%%
% parameters
%%%
Nx = 256
Lx = 2*pi;
amp_y = 1;
pulse = 2;
periodicity = 1;

%%%
% data
%%%
nx = Nx+1;
dx = Lx/(nx-1);
x_tmp = [0:dx:Lx]';
x = x_tmp(1:end-1);
y = amp_y * sin(pulse*2*pi*x/Lx+pi/2);
y_prime_th = amp_y * pulse * (2*pi/Lx) * cos(pulse*2*pi*x/Lx+pi/2);

%%%
% numerical derivatives
%%%
y_prime_df_1_p = derivative1(y,dx,1,1);
y_prime_df_2_p = derivative1(y,dx,2,1);
y_prime_df_3_p = derivative1(y,dx,3,1);
y_prime_fft = derivative1(y,dx,0,1);

%%%
% test of decomposition
%%%
delta_x = 10;

% finite differences
A(:,1)=y_prime_df_1_p;%((1+delta_x):(Nx-delta_x));
A(:,2)=y_prime_df_2_p;%((1+delta_x):(Nx-delta_x));
A(:,3)=y_prime_df_3_p;%((1+delta_x):(Nx-delta_x));
A(:,4)=y_prime_fft;%((1+delta_x):(Nx-delta_x));


% resolution of the system

w_th = eye(4);

for j = 1:10
  for i = 1:50
      
    %indicies = floor(rand(4,1)*256);%
    indicies = j + [0 i*2 i*3 i*4];
    indicies
    
    iter = 0;
    
    iter = iter +1;
    w(:,iter)= A(indicies,:)\y_prime_df_1_p(indicies);
    
    iter = iter +1;
    w(:,iter)= A(indicies,:)\y_prime_df_2_p(indicies);
    
    iter = iter +1;
    w(:,iter)= A(indicies,:)\y_prime_df_3_p(indicies);
    
    iter = iter +1;
    w(:,iter)= A(indicies,:)\y_prime_fft(indicies);
    
    w;
    
    w_err(j,i) = sum(sum(abs(w - w_th)));
    
  end
end

subplot(1,2,1)
imagesc(w_err);colorbar

subplot(1,2,2)
imagesc(log10(w_err));colorbar

if(1==0)

  reference = y_prime_fft((1+delta_x):(Nx-delta_x));
  
  subplot(2,1,1)
  plot(y_prime_df_1_p((1+delta_x):(Nx-delta_x))-reference,'r'); hold on;
  plot(y_prime_df_2_p((1+delta_x):(Nx-delta_x))-reference,'g'); hold on;
  plot(y_prime_df_3_p((1+delta_x):(Nx-delta_x))-reference,'b'); hold on;
  plot(y_prime_fft((1+delta_x):(Nx-delta_x))-reference,'k'); hold on;
  
  
  subplot(2,1,2)
  plot(log10(abs(y_prime_df_1_p((1+delta_x):(Nx-delta_x))-reference)),'r'); hold on;
  plot(log10(abs(y_prime_df_2_p((1+delta_x):(Nx-delta_x))-reference)),'g'); hold on;
  plot(log10(abs(y_prime_df_3_p((1+delta_x):(Nx-delta_x))-reference)),'b'); hold on;
  plot(log10(abs(y_prime_fft((1+delta_x):(Nx-delta_x))-reference)),'k'); hold on;
  
end
