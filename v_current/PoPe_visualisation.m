
% variables used for displaying windows
nb_run_tot = numel(data);
nb_fig_row = floor(sqrt(nb_run_tot));
nb_fig_col = ceil(nb_run_tot/nb_fig_row);

% variables used for legends
legend_300 = {};
legend_400 = {};

% display figure
figure_each_run=0;

%%%
% Analysis of PoPe results
%%%
% 1) taking into account the theoretical value as the target
for iter_run = 1:numel(data)
  for i = data(iter_run).nb_alpha
    alpha_th = data(iter_run).alpha_th;
    PoPe_c(iter_run,:) = mean(data(iter_run).alpha(:,2:end-1),2); % effective weight
    PoPe_Delta_C(iter_run,:) = abs(alpha_th - PoPe_c(iter_run,:)); % gap between effective and theoretical weight
    PoPe_delta_c(iter_run,:) = mean(abs(alpha_th'*ones(1,size(data(iter_run).alpha(:,2:end-1),2))-data(iter_run).alpha(:,2:end-1)),2);
  end
end
% 2) taking into account the last run's effective weights as the target
for iter_run = 1:numel(data)
  for i = data(iter_run).nb_alpha
    alpha_th = PoPe_c(end,:);
    PoPe_Delta_C_vs_num(iter_run,:) = abs(alpha_th - PoPe_c(iter_run,:)); % gap between effective and theoretical weight
    PoPe_delta_c_vs_num(iter_run,:) = mean(abs(alpha_th'*ones(1,size(data(iter_run).alpha(:,2:end-1),2))-data(iter_run).alpha(:,2:end-1)),2);
  end
end


for iter_run = 1:numel(data) % iteration over runs
  % iterators used to select colors
  iter_colors_1=iter_run;
  iter_colors_2=1;
  %mod(iter_run,numel(choice_derivative_range))+1;
  %floor(iter_run/numel(Nx_range))+1;
  

  %%%
  % one figure for each weigth alpha_i
  %%%

  
  for i = 1:data(iter_run).nb_alpha
    % histograms
    figure(100+i)
    subplot(nb_fig_row,nb_fig_col,iter_run)
    hist(data(iter_run).alpha(i,2:end-1)); hold on;
    title(['alpha_{i = ' num2str(i) '}, iter run = ' num2str(iter_run)])
  end

  
  %%%
  % one figure for each run
  %%%
  if(figure_each_run==1)
    figure(200+iter_run)
    subplot(3,1,1)
    imagesc(data(iter_run).time,data(iter_run).x,data(iter_run).f_save); colorbar;
    xlabel('t');
    ylabel('x');
    title(['f, iter run = ' num2str(iter_run)])
    subplot(3,1,2)
    imagesc(data(iter_run).time,data(iter_run).x,data(iter_run).dtf_save); colorbar;
    xlabel('t');
    ylabel('x');
    title(['d_t f, iter run = ' num2str(iter_run)])
    subplot(3,1,3)
    imagesc(data(iter_run).time,data(iter_run).x,data(iter_run).epsilon); colorbar;
    xlabel('t');
    ylabel('x');
    title(['epsilon, iter run = ' num2str(iter_run)])
  end

  %%%
  % one figure for all
  %%%
  
  
  figure(300);
  plot(data(iter_run).time,sum(data(iter_run).f_save),'color',colors_1(iter_colors_1,:)); hold on;
  legend_300 = { legend_300{:} ['run = ' num2str(iter_run)]}; 
  if(iter_run == numel(data))
    legend(legend_300);
    axis tight
    xlabel('t');
    ylabel('\int_x f(x) dx');
    title('conservation of f')
  end

  figure(400);
  for i = 1:data(iter_run).nb_alpha
    subplot(data(iter_run).nb_alpha,1,i)
    plot(data(iter_run).time(2:end-1),data(iter_run).alpha(i,2:end-1),'color',colors_1(iter_colors_1,:)); hold on;
    if(i==1)
      legend_400 = { legend_400{:} ['run = ' num2str(iter_run)]}; 
      if(iter_run == numel(data))
	legend_400 = { legend_400{:} 'theoretical value'}; 
      end
    end
    if(iter_run == numel(data))
      plot(data(iter_run).time(2:end-1),data(iter_run).alpha_th(i)*ones(size(data(iter_run).alpha(i,2:end-1))),'k--');
      legend(legend_400);
      xlabel('t');
      ylabel(['i = ' num2str(i)]);
      xlim([time(2) time(end-1)])      
    end
  end
   

end


for iter_run = 1:numel(data)
  scan_dx(iter_run) = data(iter_run).dx;
  scan_dt(iter_run) = data(iter_run).dt;
end




if(numel(scan_dx)>1)

%%%
% scan dx

figure(1000)
for group = 1:numel(choice_derivative_range)
    selected_runs = [group:numel(choice_derivative_range):numel(data)];
    iter_colors_2 = group;
    for i = 1:data(1).nb_alpha
      subplot(1,data(1).nb_alpha,i)
      plot(log2(scan_dx(selected_runs)),PoPe_c(selected_runs,i),'--*','color',colors_2(iter_colors_2,:)); hold on;
      title(['PoPe c, alpha_{i = ' num2str(i) '}'])
      xlabel('log2 : dx');
      ylabel('PoPe c');
      axis square; axis tight
    end
end

figure(2000)
for group = 1:numel(choice_derivative_range)
    selected_runs = [group:numel(choice_derivative_range):numel(data)];
    iter_colors_2 = group;
    for i = 1:data(1).nb_alpha
      subplot(1,data(1).nb_alpha,i)
      plot(log2(scan_dx(selected_runs)),log10(PoPe_Delta_C(selected_runs,i)),'--*','color',colors_2(iter_colors_2,:)); hold on;
      title(['PoPe Delta C, alpha_{i = ' num2str(i) '}'])
      xlabel('log2 : dx');
      ylabel('log10 : PoPe Delta C');
      axis square; axis tight
    end
end

figure(3000)
for group = 1:numel(choice_derivative_range)
    selected_runs = [group:numel(choice_derivative_range):numel(data)];
    iter_colors_2 = group;
    for i = 1:data(1).nb_alpha
      subplot(1,data(1).nb_alpha,i)
      plot(log2(scan_dx(selected_runs)),log10(PoPe_delta_c(selected_runs,i)),'--*','color',colors_2(iter_colors_2,:)); hold on;
      title(['PoPe delta c, alpha_{i = ' num2str(i) '}'])
      xlabel('log2 : dx');
      ylabel('log10 : PoPe delta c');
      axis square; axis tight
    end
end

figure(4000)
for group = 1:numel(choice_derivative_range)
    selected_runs = [group:numel(choice_derivative_range):numel(data)];
    iter_colors_2 = group;
    for i = 1:data(1).nb_alpha
      subplot(1,data(1).nb_alpha,i)
      plot(log2(scan_dx(selected_runs)),log10(PoPe_Delta_C_vs_num(selected_runs,i)),'--*','color',colors_2(iter_colors_2,:)); hold on;
      title(['PoPe Delta C, alpha_{i = ' num2str(i) '}'])
      xlabel('log2 : dx');
      ylabel('log10 : PoPe Delta C');
      axis square; axis tight
    end
end

figure(5000)
for group = 1:numel(choice_derivative_range)
    selected_runs = [group:numel(choice_derivative_range):numel(data)];
    iter_colors_2 = group;
    for i = 1:data(1).nb_alpha
      subplot(1,data(1).nb_alpha,i)
      plot(log2(scan_dx(selected_runs)),log10(PoPe_delta_c_vs_num(selected_runs,i)),'--*','color',colors_2(iter_colors_2,:)); hold on;
      title(['PoPe delta c, alpha_{i = ' num2str(i) '}'])
      xlabel('log2 : dx');
      ylabel('log10 : PoPe delta c');
      axis square; axis tight
    end
end

end % scan dx


if(numel(scan_dt)>1)

%%%
% scan dt

figure(11000)
for group = 1:numel(choice_derivative_range)
    selected_runs = [group:numel(choice_derivative_range):numel(data)];
    iter_colors_2 = group;
    for i = 1:data(1).nb_alpha
      subplot(1,data(1).nb_alpha,i)
      plot(log2(scan_dt(selected_runs)),PoPe_c(selected_runs,i),'--*','color',colors_2(iter_colors_2,:)); hold on;
      title(['PoPe c, alpha_{i = ' num2str(i) '}'])
      xlabel('log2 : dt');
      ylabel('PoPe c');
      axis square; axis tight
    end
end

figure(12000)
for group = 1:numel(choice_derivative_range)
    selected_runs = [group:numel(choice_derivative_range):numel(data)];
    iter_colors_2 = group;
    for i = 1:data(1).nb_alpha
      subplot(1,data(1).nb_alpha,i)
      plot(log2(scan_dt(selected_runs)),log10(PoPe_Delta_C(selected_runs,i)),'--*','color',colors_2(iter_colors_2,:)); hold on;
      title(['PoPe Delta C, alpha_{i = ' num2str(i) '}'])
      xlabel('log2 : dt');
      ylabel('log10 : PoPe Delta C');
      axis square; axis tight
    end
end

figure(13000)
for group = 1:numel(choice_derivative_range)
    selected_runs = [group:numel(choice_derivative_range):numel(data)];
    iter_colors_2 = group;
    for i = 1:data(1).nb_alpha
      subplot(1,data(1).nb_alpha,i)
      plot(log2(scan_dt(selected_runs)),log10(PoPe_delta_c(selected_runs,i)),'--*','color',colors_2(iter_colors_2,:)); hold on;
      title(['PoPe delta c, alpha_{i = ' num2str(i) '}'])
      xlabel('log2 : dt');
      ylabel('log10 : PoPe delta c');
      axis square; axis tight
    end
end

figure(14000)
for group = 1:numel(choice_derivative_range)
    selected_runs = [group:numel(choice_derivative_range):numel(data)];
    iter_colors_2 = group;
    for i = 1:data(1).nb_alpha
      subplot(1,data(1).nb_alpha,i)
      plot(log2(scan_dt(selected_runs)),log10(PoPe_Delta_C_vs_num(selected_runs,i)),'--*','color',colors_2(iter_colors_2,:)); hold on;
      title(['PoPe Delta C, alpha_{i = ' num2str(i) '}'])
      xlabel('log2 : dt');
      ylabel('log10 : PoPe Delta C');
      axis square; axis tight
    end
end

figure(15000)
for group = 1:numel(choice_derivative_range)
    selected_runs = [group:numel(choice_derivative_range):numel(data)];
    iter_colors_2 = group;
    for i = 1:data(1).nb_alpha
      subplot(1,data(1).nb_alpha,i)
      plot(log2(scan_dt(selected_runs)),log10(PoPe_delta_c_vs_num(selected_runs,i)),'--*','color',colors_2(iter_colors_2,:)); hold on;
      title(['PoPe delta c, alpha_{i = ' num2str(i) '}'])
      xlabel('log2 : dt');
      ylabel('log10 : PoPe delta c');
      axis square; axis tight
    end
end

end % scan dt

%   
%   figure(3000)
%   for w = 1:nb_w
%       % std
%     subplot(1,nb_w,w)
%     plot(log2(dx_run),log10(PoPe_delta_c(:,w)),'-*','color',colors_2(iter_colors_2,:));  hold on;
%     axis square; axis tight
%   end

%   
%   % %%%
%   % % scan dt
%   % 
%   % figure(1000001)
%   % for w = 1:nb_w
%   %   subplot(1,nb_w,w)
%   %   plot(log10(dt_run),PoPe_c(:,w),'-*');
%   %   title(['PoPe c, w = ' num2str(w)])
%   %   xlabel('log10 : abs(dt)');
%   %   ylabel('PoPe c');
%   %   axis square; axis tight
%   % end
%   % 
%   % figure(2000001)
%   % for w = 1:nb_w
%   %     % distance to target
%   %   subplot(1,nb_w,w)
%   %   plot(log10(dt_run),log10(PoPe_Delta_C(:,w)),'-*');
%   %   title(['PoPe Delta C, w = ' num2str(w)])
%   %   xlabel('log10 : abs(dt)');
%   %   ylabel('log10 : PoPe Delta C');
%   %   axis square; axis tight
%   % end
%   % 
%   % figure(3000001)
%   % for w = 1:nb_w
%   %     % std
%   %   subplot(1,nb_w,w)
%   %   plot(log10(dt_run),log10(PoPe_delta_c(:,w)),'-*');
%   %   title(['PoPe delta c, w = ' num2str(w)])
%   %   xlabel('log10 : abs(dt)');
%   %   ylabel('log10 : PoPe delta c');
%   %   axis square; axis tight
%   % end

%   
%   
%   
%   
%   
%   if(1==0)
%   tmp = mean(alpha(:,2:end-1),2)';
%   c_eff = tmp(1);
%   D_eff = tmp(2);
%   display(['c D : ' num2str([c D])]);
%   display(['c_eff D_eff : ' num2str([c_eff D_eff])]);
%   display(['delta_c delta_D : ' num2str([c_eff-c D_eff-D])]);
%   display(' ')
%   
%   delta_c = sqrt(mean(abs(alpha(1,2:end-1)-c).^2,2)');
%   delta_D = sqrt(mean(abs(alpha(2,2:end-1)-D).^2,2)');
%   
%   display(['delta_c delta_D : ' num2str([delta_c delta_D])]);
%   display(' ')
%   
%   sigma_c = sqrt(std(alpha(1,2:end-1)-c));
%   sigma_D = sqrt(std(alpha(2,2:end-1)-D));
%   
%   display(['sigma_c sigma_D : ' num2str([sigma_c sigma_D])]);
%   display(' ')
%   
%   
%   figure(101);
%   
%   subplot(2,2,1)
%   plot(c,c_eff,'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('c');
%   ylabel('c_{eff}');
%   plot(c_range,0*c_range,'k--');
%   plot(c_range,c_range,'k--');
%   plot(c_range,-c_range,'k--');
%   axis square
%   
%   subplot(2,2,2)
%   plot(log10(abs(c)),log10(abs(c_eff)),'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('log10 : abs(c)');
%   ylabel('log10 : abs(c_{eff})');
%   plot(log10(c_range),log10(c_range),'k--');
%   axis square
%   
%   subplot(2,2,3)
%   plot(c_eff,D_eff,'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('c_{eff}');
%   ylabel('D_{eff}');
%   axis square
%   
%   subplot(2,2,4)
%   plot(log10(abs(c_eff)),log10(abs(D_eff)),'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('log10 : abs(c_{eff})');
%   ylabel('log10 : abs(D_{eff})');
%   axis square
%   
%   figure(102);
%   
%   subplot(2,2,1)
%   plot(D,D_eff,'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('D');
%   ylabel('D_{eff}');
%   plot(D_range,0*D_range,'k--');
%   plot(D_range,D_range,'k--');
%   plot(D_range,-D_range,'k--');
%   axis square
%   
%   subplot(2,2,2)
%   plot(log10(abs(D)),log10(abs(D_eff)),'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('log10 : abs(D)');
%   ylabel('log10 : abs(D_{eff})');
%   plot(log10(D_range),log10(D_range),'k--');
%   axis square
%   
%   subplot(2,2,3)
%   plot(D_eff,c_eff,'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('D_{eff}');
%   ylabel('c_{eff}');
%   axis square
%   
%   subplot(2,2,4)
%   plot(log10(abs(D_eff)),log10(abs(c_eff)),'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('log10 : abs(D_{eff})');
%   ylabel('log10 : abs(c_{eff})');
%   axis square
%   
%   
%   figure(103)
%   
%   subplot(2,2,1)
%   plot(dx,D_eff-D,'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('dx');
%   ylabel('D_{eff}-D');
%   axis square
%   
%   subplot(2,2,2)
%   plot(log10(dx),log10(abs(D_eff-D)),'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('log10 : abs(dx)');
%   ylabel('log10 : abs(D_{eff}-D)');
%   axis square
%   
%   subplot(2,2,3)
%   plot(log10(dx),log10(delta_D),'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('log10 : abs(dx)');
%   ylabel('log10 : abs(delta_D)');
%   axis square
%   
%   subplot(2,2,4)
%   plot(log10(dx),log10(sigma_D),'*','color',colors_1(iter_colors_1,:));hold on;
%   xlabel('log10 : abs(dx)');
%   ylabel('log10 : sigma_D');
%   axis square
%   
%   end
%   
%   
%   
