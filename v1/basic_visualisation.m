
figure(1)
dt_off = 1;
dx_off = 0;
subplot(3,1,1)
imagesc(time,x,f_save(:,1+dt_off:end-dt_off)); colorbar;
xlabel('t');
ylabel('x');
title('f')
subplot(3,1,2)
imagesc(time,x,dtf_save(:,1+dt_off:end-dt_off)); colorbar
xlabel('t');
ylabel('x');
title('d_t f')
subplot(3,1,3)
dx2f_save = derivative2(f_save,dx,choice_derivative,periodicity);
D_local = dtf_save(1+dx_off:end-dx_off,1+dt_off:end-dt_off)./dx2f_save(1+dx_off:end-dx_off,1+dt_off:end-dt_off);
imagesc(time,x,D_local); colorbar
xlabel('t');
ylabel('x');
title('d_t f / d_{x2} f = D_{local}')

figure(2)
subplot(2,1,1)
plot(sum(f_save));
xlabel('t');
ylabel('\int_x f(x) dx');
title('conservation of f')
subplot(2,1,2)
hist(D_local(:));
title('hist D_{local}')

figure(3)
fft_f_save = fft(f_save);
subplot(2,1,1)
imagesc(time,[0:Nx/2-1],abs(fft_f_save(1:Nx/2,:)));
xlabel('t');
ylabel('mode');
title('abs fft f')
subplot(2,1,2)
plot(time,abs(fft_f_save(1+pulse,:)));
xlabel('t');
ylabel(['f_{fft}(n=' num2str(pulse) ')']);

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
