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
delta_x = 0;

w_th = eye(4);

for jump_x = 1:50;

clear A

% finite differences
A(:,1)=y_prime_df_1_p((1+delta_x):jump_x:(Nx-delta_x));
A(:,2)=y_prime_df_2_p((1+delta_x):jump_x:(Nx-delta_x));
A(:,3)=y_prime_df_3_p((1+delta_x):jump_x:(Nx-delta_x));
A(:,4)=y_prime_fft((1+delta_x):jump_x:(Nx-delta_x));

A_mc = A'*A;

% resolution of the system

iter = 0;

iter = iter +1;
w(:,iter)= A_mc\(A'*y_prime_df_1_p((1+delta_x):jump_x:(Nx-delta_x)));

iter = iter +1;
w(:,iter)= A_mc\(A'*y_prime_df_2_p((1+delta_x):jump_x:(Nx-delta_x)));

iter = iter +1;
w(:,iter)= A_mc\(A'*y_prime_df_3_p((1+delta_x):jump_x:(Nx-delta_x)));

iter = iter +1;
w(:,iter)= A_mc\(A'*y_prime_fft((1+delta_x):jump_x:(Nx-delta_x)));

  w_err(jump_x) = sum(sum(abs(w - w_th)));
    
end

subplot(2,2,1)
plot(w_err)
subplot(2,2,3)
plot(log10(w_err))

[p v] = hist(w_err,10)

subplot(2,2,2)
plot(v,p)

[p v] = hist(log10(w_err),10)

subplot(2,2,4)
plot(v,p)

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
