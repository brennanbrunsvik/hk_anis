clear; 
clc; 
addpath('~/Documents/repositories/Base_code/colormaps/colormaps_ander_biguri'); 

%% Prepare model. 
load('models/model_ta_kmsc.mat'); 
sta_name = 'TA.KMSC'; 

% Sort by ray parameter. Makes plotting more obvious. 
[rayp, sortp] = sort(rayp); 
rf = rf(:,sortp); 

phi = 1./xi; % Flip, get h on top and v on bottom

% Get average crustal values to be consistent with Zhu and Kanamori math
[rhoav, vsav, vpav, xiav, phiav, etaav] =  hk_average_crust(...
    rho, vs, vp, xi, phi, eta, z, zmoh); 

%% Loop over receiver functions making HK stack for each. 
[Exi, H, K, t_pred_all_xi] = hk_anis_loop(...
    rf, tt, rayp, vsav, rhoav, xiav, phiav, etaav); % accounting for xi
[E00, H, K, t_pred_all_00] = hk_anis_loop(...
    rf, tt, rayp, vsav, rhoav, 1, 1, 1); % not accounting for anisotropy

%% Figure out relevant predicted pulse timings. 
[ihmax, ikmax] = find(Exi == max(Exi ,[], 'all')); 
kmax = K(ikmax); 
hmax = H(ihmax); 

t_pred_xi = hk_pulse_time_interp(H, K, t_pred_all_xi, hmax, kmax); 
t_pred_00 = hk_pulse_time_interp(H, K, t_pred_all_00, hmax, kmax); 


%% Plots. 
% HK stack 
fhand_norm = @(inval)inval ./ max(max(inval)); % Return normalized inval 

figure(301); clf; hold on; set(gcf, 'pos', [-1089 329 364 218]); 
subplot(1,1,1); hold on; 
set(gca,'ydir', 'reverse', 'LineWidth', 1.5);
grid on; 
box on; 
xlabel('\kappa'); 
ylabel('H (km)'); 
title(sprintf('H-\\kappa stack %s', sta_name), 'fontweight', 'normal'); 
% contourf(K, H, Exi_all{i_xi_true}', 30, 'EdgeAlpha', 0.1); 

plt_ylim = [27, 40]; 
plt_xlim = [1.6, 2.1]; 
ylim(plt_ylim); 
xlim(plt_xlim); 

lvl_cnt = [0.7, 0.95]; 
LW = 1; 
[~,hnd_xistart] = contour(K, H,  Exi,...
    lvl_cnt, 'k', 'LineWidth', LW*1.5, ...
    'DisplayName', sprintf('\\xi = %1.2f', xiav ) ); 
[~,hnd_xiend  ] = contour(K, H,  E00,...
    lvl_cnt, 'r', 'LineWidth', LW*1.5, ...
    'DisplayName', sprintf('\\xi = %1.2f', 1.00 )); 

legend([hnd_xistart, hnd_xiend], 'Location', 'best'); 

exportgraphics(gcf, sprintf('figs/hk_with_without_anis_%s.pdf',sta_name), ...
    'ContentType', 'vector'); 

%% Nice looking HK stack
figure(302); clf; hold on; set(gcf, 'pos', [2623 611 325 212]); 
subplot(1,1,1); hold on; 
% ax1 = gca(); 
% ax2 = copyobj(gca, gcf); 
set(gca,'ydir', 'reverse', 'LineWidth', 1.5);
box on; 


xlabel('\kappa'); 
ylabel('H (km)'); 
title(sprintf('H-\\kappa stack %s', sta_name), 'fontweight', 'normal'); 


[~,contH] = contourf(gca, K, H, Exi, 15, 'linestyle', 'none'); 

cbar = colorbar(gca,'eastoutside'); 
set(cbar, 'fontsize', 12); 

try 
    colormap(viridis()); 
catch 
    warning('Missing colormap viridis. Should be in repositories somewhere. '); 
end

scatter(kmax, hmax, 150, 'r', 'pentagram'); 
% linkaxes([ax1, ax2]); 

exportgraphics(gcf, sprintf('figs/hk_nice_%s.jpeg',sta_name),...
    'Resolution', 700); 

%% Receiver function
figure(202); clf; hold on; 
set(gcf, 'pos', [2032 336 578 1029]); 
set(gca, 'LineWidth', 1.5, 'XGrid', 'on', 'XMinorTick', 'on'); box on; %grid on; 
xlabel('Time (s)'); 
title('Phase timing', 'FontWeight','normal'); 
set(gca, 'YTick', []); 
xlim([-1, 20])
yshift_const = 0.05; 

nrf = size(rf,2); 

for irf = 1:nrf
    yshift = irf * yshift_const; 
    rfi = rf(:,irf);

    t_plot_xi = t_pred_xi(irf,:); 
    t_plot_00 = t_pred_00(irf,:)'; 

    hnd_t_xi = scatter(...
        t_plot_xi', yshift + interp1(tt, rf(:,irf), t_plot_xi, 'linear'),...
        50, 'blue', '+', 'LineWidth', 2, 'DisplayName', 'With \xi'); % If using true parameters and anisotropic stack
    hnd_t_00 = scatter(...
        t_plot_00', yshift + interp1(tt, rf(:,irf), t_plot_00, 'linear'),...
        50, 'red', '+', 'LineWidth', 2, 'DisplayName', 'No \xi'); % If using true parameters and isotropic stack
    hnd_rf = plot(tt, yshift+rfi, 'k', 'linewidth', 0.75);
end

ylim([-10*yshift_const, yshift_const * (nrf+10)])
lgd = legend([hnd_t_00, hnd_t_xi], 'Location', 'southeast'); 

exportgraphics(gcf, sprintf('figs/rfs_%s.pdf',sta_name), ...
    'ContentType', 'vector'); 