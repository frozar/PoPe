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
y_prime_df_1_nop = derivative1(y,dx,1,0);
y_prime_df_1_p = derivative1(y,dx,1,1);
y_prime_df_2_nop = derivative1(y,dx,2,0);
y_prime_df_2_p = derivative1(y,dx,2,1);
y_prime_df_3_nop = derivative1(y,dx,3,0);
y_prime_df_3_p = derivative1(y,dx,3,1);
y_prime_fft = derivative1(y,dx,0,1);
y_sec_fft = derivative2(y,dx,0,1);


%%%
% test of decomposition
%%%
delta_x = 10;
iter = 0;

% finite differences
if(1==0)
  % non periodic
  iter = iter+1;
  A(:,iter)=y_prime_df_1_nop((1+delta_x):(Nx-delta_x));
  iter = iter+1;
  A(:,iter)=y_prime_df_2_nop((1+delta_x):(Nx-delta_x));
  iter = iter+1;
  A(:,iter)=y_prime_df_3_nop((1+delta_x):(Nx-delta_x));
else
  % periodic
  iter = iter+1;
  A(:,iter)=y_prime_df_1_p((1+delta_x):(Nx-delta_x));
  iter = iter+1;
  A(:,iter)=y_prime_df_2_p((1+delta_x):(Nx-delta_x));
  iter = iter+1;
  A(:,iter)=y_prime_df_3_p((1+delta_x):(Nx-delta_x));
end

if(1==0)
  % fft
  iter = iter+1;
  A(:,iter)=y_prime_fft((1+delta_x):(Nx-delta_x));
end

if(1==0)
  % second derivative fft 
  iter = iter+1;
  A(:,iter)=y_sec_fft((1+delta_x):(Nx-delta_x));
end

A1(:,1) = y_prime_df_1_p((1+delta_x):(Nx-delta_x));
A1(:,2) = 1;

A2(:,1) = y_prime_df_2_p((1+delta_x):(Nx-delta_x));
A2(:,2) = 1;

A3(:,1) = y_prime_df_3_p((1+delta_x):(Nx-delta_x));
A3(:,2) = 1;

iter = 0;

if(1==0)
  
  iter = iter +1;
  w(:,iter)=A\y_prime_df_1_nop((1+delta_x):(Nx-delta_x));
  
  iter = iter +1;
  w(:,iter) = A\y_prime_df_2_nop((1+delta_x):(Nx-delta_x));
  
  iter = iter +1;
  w(:,iter) = A\y_prime_df_3_nop((1+delta_x):(Nx-delta_x));

else
  
  iter = iter +1;
  w(:,iter)= A\y_prime_df_1_p((1+delta_x):(Nx-delta_x));
  
  iter = iter +1;
  w(:,iter)= A\y_prime_df_2_p((1+delta_x):(Nx-delta_x));
  
  iter = iter +1;
  w(:,iter)= A\y_prime_df_3_p((1+delta_x):(Nx-delta_x));

end

w

w1(:,1)= A1\y_prime_df_1_p((1+delta_x):(Nx-delta_x));
w1(:,2)= A2\y_prime_df_1_p((1+delta_x):(Nx-delta_x));
w1(:,3)= A3\y_prime_df_1_p((1+delta_x):(Nx-delta_x));

w1

w2(:,1)= A1\y_prime_df_2_p((1+delta_x):(Nx-delta_x));
w2(:,2)= A2\y_prime_df_2_p((1+delta_x):(Nx-delta_x));
w2(:,3)= A3\y_prime_df_2_p((1+delta_x):(Nx-delta_x));

w2

w3(:,1)= A1\y_prime_df_3_p((1+delta_x):(Nx-delta_x));
w3(:,2)= A2\y_prime_df_3_p((1+delta_x):(Nx-delta_x));
w3(:,3)= A3\y_prime_df_3_p((1+delta_x):(Nx-delta_x));

w3
