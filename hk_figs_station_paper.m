clear; 
clc; 
restoredefaultpath; 

addpath("elastic/"); 
addpath("hk_calculation/"); 
addpath("hk_calculation/"); 
addpath('~/Documents/repositories/Base_code/colormaps/colormaps_ander_biguri'); 

%% Prepare model. 
load('models/model_ta_kmsc.mat'); 
sta_name = 'TA.KMSC'; 

%% Remove bad receiver functions

% First row: Handle bad parent pulse. Second is for whole receiver function. 
% First collumn: Cutoff for noramlized average cross correlation to other
% receiver functions. Second collumn: Time window (absolute value)
acparms = [0.95, 1; 0.55, inf]; % auto correlation parameters. These don't need to be perfect. It's just removing really bad receiver functions. 

for iacremove = 1:2; 
    nrf = size(rf,2); % Changes with each iacremove
    ac = nan(nrf, nrf); % Autocorrelation matrix. 
    twin = abs(tt) < acparms(iacremove,2); % Window within which is the parent pulse, to autocorrelate. 
    ac_cutoff = acparms(iacremove,1); 
    for irf = 1:nrf; 
        for jrf = 1:nrf; % This is fast and done not often so don't worry about taking advantage of symmetry
            rfi = rf(twin,irf); 
            rfj = rf(twin,jrf); 
            ac(irf, jrf) = rfi' * rfj ./ ...
                ( sqrt(rfi' * rfi)*sqrt(rfj' * rfj) ); 
        end
    end
    acs = sum(ac) ./ length(ac); % autocorrelation sum. Normalize it. 
    ac_keep = acs > ac_cutoff; 
    
    rayp     = rayp(:,ac_keep     ); 
    rf       = rf  (:,ac_keep     ); 
end

%% Sort by ray parameter. Makes plotting more obvious. 
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
    rf, tt, rayp, vsav, rhoav,    1,     1,     1); % not accounting for anisotropy

%% Figure out relevant predicted pulse timings. 
[ihmax, ikmax] = find(Exi == max(Exi ,[], 'all')); % If you were to optimize while accounting for anisotropy
kmax = K(ikmax); 
hmax = H(ihmax); 

[ihmaxiso, ikmaxiso] = find(E00 == max(E00 ,[], 'all')); % If you were to optimize while ignoring anisotropy
kmaxiso = K(ikmaxiso); 
hmaxiso = H(ihmaxiso); 

t_pred_xi = hk_pulse_time_interp(H, K, t_pred_all_xi, hmax, kmax); % using xi
t_pred_00 = hk_pulse_time_interp(H, K, t_pred_all_00, hmax, kmax); % O anisotropy, 
t_pred_0b = hk_pulse_time_interp(H, K, t_pred_all_00, hmaxiso, kmaxiso); % 0 anisotropy but best model

%% Set figure layout
c_xilow = [159, 16, 199]./255; 
c_xihi  = [17, 130, 32]./255; 
c_mod_true = 'y'; % Yellow
c_t_with_xi = [3, 15, 252]./255; ;% Anisotropic
c_t_no_xi = [255, 0, 0]./255; % not a solution at all. xi=1 but h_xi and k_xi
c_t_iso = [0, 255, 0]./255; % What would be got if people use isotropic stack

figmain = figure(502); clf; 
set(figmain, 'pos', [2476 754 810 245]); 
nxt = 5; 
nyt = 10;
xy_to_t = @(x,y)(y-1)*nyt+x; 
tiledlayout(nxt, nyt, 'TileSpacing','compact');
ax1 = nexttile(xy_to_t(1,1), [5,6]); hold on; 
ax2 = nexttile(xy_to_t(7,1), [5,4]); hold on; 

each_ax = [ax1, ax2]; 
for ieach_ax = 1:length(each_ax); 
    axes(each_ax(ieach_ax)); 
    set(gca, 'LineWidth', 1.5); box on; 
end

axes(ax1); 
xlabel('Time (s)'); 
yticklabels([]); 
title('Phase timing', 'FontWeight','normal'); 

axes(ax2); 
xlabel('\kappa'); 
ylabel('H (km)'); % yticklabels([]); 
grid on;
title(sprintf('HK stack %s', sta_name), 'fontweight', 'normal'); 


% HK stack 
axes(ax2); 
fhand_norm = @(inval)inval ./ max(max(inval)); % Return normalized inval 
set(gca,'ydir', 'reverse');

plt_ylim = [29, 39]; 
plt_xlim = [1.6, 2.1]; 
ylim(plt_ylim); 
xlim(plt_xlim); 

lvl_cnt = [0.7, 0.95]; 
LW = 1; 

LW_scale = 1.75; 
[~,hnd_xistart] = contour(K, H,  Exi,...
    lvl_cnt, 'LineWidth', LW*LW_scale, 'color', c_xihi, ...
    'DisplayName', sprintf('\\xi = %1.2f', xiav ) ); 
[~,hnd_xiend  ] = contour(K, H,  E00,...
    lvl_cnt, 'k', 'LineWidth', LW*LW_scale, ...
    'DisplayName', sprintf('\\xi = %1.2f', 1.00 )); 

% hnd_mxi  = scatter(kmax   , hmax   , 100, 'd', 'filled', ...
%     'LineWidth', 1, 'MarkerEdgeColor', 'k',... 
%     'DisplayName', 'H_{\xi},\kappa_{\xi}', 'MarkerFaceColor', c_t_with_xi ); 
% hnd_miso = scatter(kmaxiso, hmaxiso, 100, 'd', 'filled', ...
%     'LineWidth', 1, 'MarkerEdgeColor', 'k',... 
%     'DisplayName', 'H_1,\kappa_{1}', 'MarkerFaceColor', c_t_iso); 

% legend([hnd_xistart, hnd_xiend, hnd_mxi, hnd_miso], 'Location', 'best'); 
legend([hnd_xistart, hnd_xiend], 'Location', 'best'); 


% exportgraphics(gcf, sprintf('figs/hk_with_without_anis_%s.pdf',sta_name), ...
%     'ContentType', 'vector'); 

% Receiver function no ysfhit
axes(ax1); 
set(gca, 'XGrid', 'on', 'XMinorTick', 'on'); 
xlim([-1, 20])
yshift_const = 0.00; 

txt_t_xi = '$t(H_{\xi}, \kappa_{\xi}, \xi)$'; 
txt_t_0b = '$t(H_{1}, \kappa_{1}, 1)$'; 
txt_t_00 = '$t(H_{\xi}, \kappa_{\xi}, 1)$'; 

nrf = size(rf,2); 

for irf = 1:nrf
    yshift = irf * yshift_const; 
    rfi = rf(:,irf);
    hnd_rf = plot(tt, yshift+rfi, 'linewidth', 0.75,...
        'Color', [0,0,0,0.14]); ; 
end

for irf = 1:nrf
    yshift = irf * yshift_const; 
    rfi = rf(:,irf);

    t_plot_xi = t_pred_xi(irf,:); 
    t_plot_00 = t_pred_00(irf,:)'; 
    t_plot_0b = t_pred_0b(irf,:)'; 

    LWS = 1; 
    sct_symb = '.'; 
    hnd_t_xi = scatter(...
        t_plot_xi', yshift + interp1(tt, rf(:,irf), t_plot_xi, 'linear'),...
        50, sct_symb, 'LineWidth', LWS, 'DisplayName', txt_t_xi,...
        'CData', c_t_with_xi); % If using true parameters and anisotropic stack
    hnd_t_00 = scatter(...
        t_plot_00', yshift + interp1(tt, rf(:,irf), t_plot_00, 'linear'),...
        50, sct_symb, 'LineWidth', LWS, 'DisplayName', txt_t_00,...
        'CData', c_t_no_xi); % If using true parameters and isotropic stack
    hnd_t_0b = scatter(...
        t_plot_0b', yshift + interp1(tt, rf(:,irf), t_plot_0b, 'linear'),...
        50, sct_symb, 'LineWidth', LWS, 'DisplayName', txt_t_0b,...
        'CData', c_t_iso); % If using parameters from optimizing the isotropic stack
end

% hnd_t_xi = copyobj(hnd_t_xi, gca); hnd_t_xi.SizeData = 1000; 
% hnd_t_00 = copyobj(hnd_t_00, gca); hnd_t_00.SizeData = 1000; 
% hnd_t_0b = copyobj(hnd_t_0b, gca); hnd_t_0b.SizeData = 1000; 


ylim([-0.25, 0.35]); 
lgd = legend([hnd_t_xi, hnd_t_0b, hnd_t_00], 'Location', 'northeast',...
    'Interpreter','latex', 'fontsize', 12); 

exportgraphics(gcf, sprintf('figs/rfs_paper_merged_%s.pdf',sta_name), ...
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
    t_plot_0b = t_pred_0b(irf,:)'; 

    hnd_t_xi = scatter(...
        t_plot_xi', yshift + interp1(tt, rf(:,irf), t_plot_xi, 'linear'),...
        50, '+', 'CData', c_t_with_xi, 'LineWidth', 2, 'DisplayName', txt_t_xi); % If using true parameters and anisotropic stack
    hnd_t_00 = scatter(...
        t_plot_00', yshift + interp1(tt, rf(:,irf), t_plot_00, 'linear'),...
        50, '+', 'CData', c_t_no_xi, 'LineWidth', 2, 'DisplayName', txt_t_00); % If using true parameters and isotropic stack
    hnd_t_0b = scatter(...
        t_plot_0b', yshift + interp1(tt, rf(:,irf), t_plot_0b, 'linear'),...
        50, '+', 'CData', c_t_iso, 'LineWidth', 2, 'DisplayName', txt_t_0b); % If using parameters from optimizing the isotropic stack

    hnd_rf = plot(tt, yshift+rfi, 'k', 'linewidth', 0.75);
end

ylim([-10*yshift_const, yshift_const * (nrf+10)])
lgd = legend([hnd_t_xi, hnd_t_0b, hnd_t_00], 'Location', 'southeast',...
    'Interpreter','latex', 'fontsize', 13); 

text(1, -3*yshift_const     , sprintf('$p=%1.2f$',rayp(1  )), 'interpreter', 'latex' ); % State ray parameter. If it's seconds per degree, should have something like 4.5 to 9. Seconds per kilometer should be something like 0.04 to 0.08 (See Zhu and Kanamori 2000). 
text(1, (nrf+5)*yshift_const, sprintf('$p=%1.2f$',rayp(end)), 'interpreter', 'latex'  ); 

exportgraphics(gcf, sprintf('figs/rfs_%s.pdf',sta_name), ...
    'ContentType', 'vector'); 

% % % %% Receiver function no ysfhit
% % % figure(203); clf; hold on; 
% % % set(gcf, 'pos', [2032 336 578 229]); 
% % % set(gca, 'LineWidth', 1.5, 'XGrid', 'on', 'XMinorTick', 'on'); box on; %grid on; 
% % % xlabel('Time (s)'); 
% % % title('Phase timing', 'FontWeight','normal'); 
% % % set(gca, 'YTick', []); 
% % % xlim([-1, 20])
% % % yshift_const = 0.00; 
% % % 
% % % nrf = size(rf,2); 
% % % 
% % % for irf = 1:nrf
% % %     yshift = irf * yshift_const; 
% % %     rfi = rf(:,irf);
% % % 
% % %     t_plot_xi = t_pred_xi(irf,:); 
% % %     t_plot_00 = t_pred_00(irf,:)'; 
% % %     t_plot_0b = t_pred_0b(irf,:)'; 
% % % 
% % %     hnd_t_xi = scatter(...
% % %         t_plot_xi', yshift + interp1(tt, rf(:,irf), t_plot_xi, 'linear'),...
% % %         50, 'blue', '+', 'LineWidth', 2, 'DisplayName', txt_t_xi); % If using true parameters and anisotropic stack
% % %     hnd_t_00 = scatter(...
% % %         t_plot_00', yshift + interp1(tt, rf(:,irf), t_plot_00, 'linear'),...
% % %         50, 'green', '+', 'LineWidth', 2, 'DisplayName', txt_t_00); % If using true parameters and isotropic stack
% % %     hnd_t_0b = scatter(...
% % %         t_plot_0b', yshift + interp1(tt, rf(:,irf), t_plot_0b, 'linear'),...
% % %         50, 'red', '+', 'LineWidth', 2, 'DisplayName', txt_t_0b); % If using parameters from optimizing the isotropic stack
% % %     hnd_rf = plot(tt, yshift+rfi, 'k', 'linewidth', 0.75);
% % % end
% % % 
% % % 
% % % ylim([-0.25, 0.5])
% % % lgd = legend([hnd_t_xi, hnd_t_0b, hnd_t_00], 'Location', 'northeast',...
% % %     'Interpreter','latex', 'fontsize', 12); 
% % % 
% % % exportgraphics(gcf, sprintf('figs/rfs_noyshift_%s.pdf',sta_name), ...
% % %     'ContentType', 'vector'); 