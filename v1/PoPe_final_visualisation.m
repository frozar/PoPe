nb_w = size(alpha,1);

figure(1000000)
for w = 1:nb_w
  subplot(1,nb_w,w)
  plot(log10(dx_run),PoPe_c(:,w));
  title(['PoPe c, w = ' num2str(w)])
  axis square; axis tight
end

figure(400000)
for w = 1:nb_w
    % distance to target
  subplot(1,nb_w,w)
  plot(log10(dx_run),log10(PoPe_Delta_C(:,w)));
  title(['PoPe Delta C, w = ' num2str(w)])
  axis square; axis tight
end

figure(500000)
for w = 1:nb_w
    % std
  subplot(1,nb_w,w)
  plot(log10(dx_run),log10(PoPe_delta_c(:,w)));
  title(['PoPe delta c, w = ' num2str(w)])
  axis square; axis tight
end

figure(1000001)
for w = 1:nb_w
  subplot(1,nb_w,w)
  plot(log10(dt_run),PoPe_c(:,w));
  title(['PoPe c, w = ' num2str(w)])
  axis square; axis tight
end

figure(400001)
for w = 1:nb_w
    % distance to target
  subplot(1,nb_w,w)
  plot(log10(dt_run),log10(PoPe_Delta_C(:,w)));
  title(['PoPe Delta C, w = ' num2str(w)])
  axis square; axis tight
end

figure(500001)
for w = 1:nb_w
    % std
  subplot(1,nb_w,w)
  plot(log10(dt_run),log10(PoPe_delta_c(:,w)));
  title(['PoPe delta c, w = ' num2str(w)])
  axis square; axis tight
end
