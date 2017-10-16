%function PoPe_visualisation(alpha,iter_run,time,x,f_save,dtf_save,epsilon,colors,w_th)

nb_w = size(alpha,1);
iter_color = iter_run;

%%%
% one figure for each weigth

for w = 1:nb_w
  % histograms
  figure(1000+w)
  subplot(2,2,iter_run)
  hist(alpha(w,2:end-1)); hold on;
  title(['w = ' num2str(w) ', iter run = ' num2str(iter_run)])
end

%%%
% one figure for each run

figure(10000+iter_run)
subplot(2,1,1)
imagesc(time,x,dtf_save); colorbar;
xlabel('t');
ylabel('x');
title(['d_t f, iter run = ' num2str(iter_run)])
subplot(2,1,2)
imagesc(time,x,epsilon); colorbar;
xlabel('t');
ylabel('x');
title(['epsilon, iter run = ' num2str(iter_run)])

%%%
% one figure for all

figure(100000)
plot(time,sum(f_save),'color',colors(iter_color,:)); hold on;
xlabel('t');
ylabel('\int_x f(x) dx');
title('conservation of f')
xlim([time(1) time(end)])

figure(200000);
for w = 1:nb_w
  subplot(nb_w,1,w)
  plot(time(2:end-1),alpha(w,2:end-1),'color',colors(iter_color,:)); hold on;
  plot(time(2:end-1),w_th(w)*ones(size(alpha(w,2:end-1))),'k-*');
  xlabel('t');
  ylabel(['w = ' num2str(w)]);
  xlim([time(2) time(end-1)])
end

if(1==0)
tmp = mean(alpha(:,2:end-1),2)';
c_eff = tmp(1);
D_eff = tmp(2);
display(['c D : ' num2str([c D])]);
display(['c_eff D_eff : ' num2str([c_eff D_eff])]);
display(['delta_c delta_D : ' num2str([c_eff-c D_eff-D])]);
display(' ')

delta_c = sqrt(mean(abs(alpha(1,2:end-1)-c).^2,2)');
delta_D = sqrt(mean(abs(alpha(2,2:end-1)-D).^2,2)');

display(['delta_c delta_D : ' num2str([delta_c delta_D])]);
display(' ')

sigma_c = sqrt(std(alpha(1,2:end-1)-c));
sigma_D = sqrt(std(alpha(2,2:end-1)-D));

display(['sigma_c sigma_D : ' num2str([sigma_c sigma_D])]);
display(' ')


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


figure(103)

subplot(2,2,1)
plot(dx,D_eff-D,'*','color',colors(iter_color,:));hold on;
xlabel('dx');
ylabel('D_{eff}-D');
axis square

subplot(2,2,2)
plot(log10(dx),log10(abs(D_eff-D)),'*','color',colors(iter_color,:));hold on;
xlabel('log10 : abs(dx)');
ylabel('log10 : abs(D_{eff}-D)');
axis square

subplot(2,2,3)
plot(log10(dx),log10(delta_D),'*','color',colors(iter_color,:));hold on;
xlabel('log10 : abs(dx)');
ylabel('log10 : abs(delta_D)');
axis square

subplot(2,2,4)
plot(log10(dx),log10(sigma_D),'*','color',colors(iter_color,:));hold on;
xlabel('log10 : abs(dx)');
ylabel('log10 : sigma_D');
axis square

end
