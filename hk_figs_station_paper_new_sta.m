% Do anisotropic HK stacks for a particular station. 
% Make figures for use in the paper. 

% TODO I used the phrase "autocorrelation" but mean zero-lag
% cross-correlation
% TODO plotting receiver function pulse timings is slow. 

clear; 
clc; 
restoredefaultpath; 

addpath("elastic/"); 
addpath("hk_calculation/"); 
addpath("hk_calculation/"); 
addpath('~/Documents/repositories/Base_code/colormaps/colormaps_ander_biguri'); 

%% Prepare model. 
sta_name = 'AE.X16A'; 
f = load("models/download_rfs/Ears/gauss_2.5/"+sta_name+"/rfArr.mat"); 
rf = f.rf; tt = f.tt; rayParmSecDeg = f.rayParmSecDeg; incAngP = f.incAngP; 
lat = f.stla; lon = f.stlo; 
rayp = rayParmSecDeg; 
% load('models/model_ta_kmsc.mat'); 

xi = 1.33; 
VP = 6.276; %LB BMN values, mostly the same as for US ELK
VS = 1/1.76 * VP; % Specifying because I don't remember which we hold constant. 
zmoh = 34; % ELK estimated at 30 km. Use this for some rough approximations, not so important. 
rhoest = 2.83; 
% eta = 1; 
eta = 1; 
phi_xi_const = 1; % define 1/phi = phi_xi_const * xi. Still should double check how this is implimented in code below, brb2023.03.21

kNum = 510; 
hNum = 500; 
% xiNum = 31; 
plot_true = true; 
hBounds = [10, 55]; 

%% Remove bad receiver functions. 
% Remove those that have small zero-lag cross-correlation to other receiver functions. 

% auto correlation parameters. These don't need to be perfect. It's just removing really bad receiver functions. 
% First row: Handle bad parent pulse. Second is for whole receiver function. 
% First collumn: Cutoff for noramlized average cross correlation to other
% receiver functions. Second collumn: Time window (absolute value)
acparms = [0.95, 1; 0.55, inf]; 

for iacremove = 1:2; % Loop over checking the first few seconds of the receiver function, and the whole receiver function. 
    nrf = size(rf,2); % Number of receiver functions. Changes with each iacremove
    ac = nan(nrf, nrf); % Autocorrelation matrix. 
    twin = abs(tt) < acparms(iacremove,2); % Window within which is the parent pulse, to cross-correlate.  
    ac_cutoff = acparms(iacremove,1); % Get parameters for this iacremove. 
    for irf = 1:nrf; % Loop over each receiver function
        for jrf = 1:nrf; % Again loop over each receiver function. This is fast enough and done not often so don't worry about taking advantage of symmetry
            rfi = rf(twin,irf); 
            rfj = rf(twin,jrf); 
            ac(irf, jrf) = rfi' * rfj ./ ...
                ( sqrt(rfi' * rfi)*sqrt(rfj' * rfj) ); % cross-correlation. 
        end
    end
    acs = sum(ac) ./ length(ac); % Cross-correlation sum. Normalize it. 
    ac_keep = acs > ac_cutoff; % Only keep if there is high-enough average cross-correlation. 
    
    rayp     = rayp(:,ac_keep     ); % Remove ray parameter for receiver functions we aren't keeping. 
    rf       = rf  (:,ac_keep     ); % Remove bad receiver functions. 
end

%% Sort by ray parameter. Makes plotting more obvious. 
[rayp, sortp] = sort(rayp); 
rf = rf(:,sortp); 

z = 1:100; 
z = z'; 
xi = ones(size(z)) .* xi; 
eta = ones(size(z)) .* eta; 
phi = 1./phi_xi_const .* 1./xi; % See paper for description of this relationship. 
nz = length(z); 
sz = size(z); 
rho = rhoest * ones(sz); 
vs = VS * ones(sz); 
vp = VP * ones(sz); 

% Get average crustal values to use with with Zhu and Kanamori math
[rhoav, vsav, vpav, xiav, phiav, etaav] =  hk_average_crust(...
    rho, vs, vp, xi, phi, eta, z, zmoh); 

%% Loop over receiver functions making HK stack for each. 
% Note: Exi, include xi. E00, standard process without xi.  

[Exi, H, K, t_pred_all_xi] = hk_anis_loop(...
    rf, tt, rayp, vsav, rhoav, xiav, phiav, etaav, ...
            'kNum', kNum, 'hNum', hNum, 'normalize', false, ...
            'hBounds', hBounds); % accounting for xi
[E00, H, K, t_pred_all_00] = hk_anis_loop(...
    rf, tt, rayp, vsav, rhoav,    1,     1,     1, ...
            'kNum', kNum, 'hNum', hNum, 'normalize', false, ...
            'hBounds', hBounds); ; % not accounting for anisotropy

%% Figure out relevant predicted pulse timings. 
[ihmax, ikmax] = find(Exi == max(Exi ,[], 'all')); % If you were to optimize while accounting for anisotropy, which k and h have the highest energy?
kmax = K(ikmax); 
hmax = H(ihmax); 

[ihmaxiso, ikmaxiso] = find(E00 == max(E00 ,[], 'all')); % If you were to optimize while ignoring anisotropy
kmaxiso = K(ikmaxiso); 
hmaxiso = H(ihmaxiso); 

Exi_max = max(Exi, [], 'all'); 
E00_max = max(E00, [], 'all'); 


t_pred_xi = hk_pulse_time_interp(H, K, t_pred_all_xi, hmax, kmax); % using xi
t_pred_00 = hk_pulse_time_interp(H, K, t_pred_all_00, hmax, kmax); % O anisotropy, 
t_pred_0b = hk_pulse_time_interp(H, K, t_pred_all_00, hmaxiso, kmaxiso); % 0 anisotropy but best model

%% Start figures. 

% Some custom colors. 
c_xilow = [179, 54, 214]./255; % Low values of xi
c_xihi = [201, 135, 2]./255; 
c_mod_true = 'y'; % Yellow, star? 
c_t_with_xi = [3, 15, 252]./255; ;% Anisotropic 
c_t_no_xi = [230, 2, 14]./255; 
c_t_iso = [10, 247, 30]./255; 
cnt_mod = 1/1.5; % Multiply colors by this when doing contours. Make them a bit darker. 


figmain = figure(502); clf; 
set(figmain, 'pos', [2476 754 810 245]); 
nxt = 5; % x tiles for tiledlayout
nyt = 10; % y tiles for tiledlayout
xy_to_t = @(x,y)(y-1)*nyt+x; % A hack so we can fill multiple tiles with one figure. 
tiledlayout(nxt, nyt, 'TileSpacing','compact');
ax1 = nexttile(xy_to_t(1,1), [5,6]); hold on; 
ax2 = nexttile(xy_to_t(7,1), [5,4]); hold on; 

each_ax = [ax1, ax2]; % Prettify
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

lvl_cnt = [0.8, 0.95] ; % .* max([Exi, E00], [], 'all'); % Contour percentages for HK stack plot. Level is percent compared to maximum of either HK stack. 
LW = 1; % Linewidth. 

LW_scale = 1.75; % Increase linewidth by this much for contour lines (?) TODO verify
[~,hnd_xistart] = contour(K, H,  Exi,...
    lvl_cnt * max(Exi, [], 'all'), 'LineWidth', LW*LW_scale, 'color', c_t_with_xi.*cnt_mod, ...
    'DisplayName', sprintf('\\xi = %1.2f', xiav ) ); 
[~,hnd_xiend  ] = contour(K, H,  E00,...
    lvl_cnt * max(E00, [], 'all'), 'LineWidth', LW*LW_scale, 'color', c_t_iso.*cnt_mod, ...
    'DisplayName', sprintf('\\xi = %1.2f', 1.00 )); 
% text(kmax   , hmax   , sprintf('E=\n%1.2f',Exi_max/E00_max), 'color', c_xihi, ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); 
% text(kmaxiso, hmaxiso, sprintf('E=\n1.00'), 'color', 'k', ...
%     'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle'); 
text(0.02, 0.98, ['E_{\xi} / E_{1} = ' sprintf('%1.2f', Exi_max/E00_max)],...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Units', 'normalized'); 
scatter(kmax, hmax, 20, c_t_with_xi.*cnt_mod, 'filled'); 
scatter(kmaxiso, hmaxiso, 20, c_t_iso.*cnt_mod, 'filled'); 

legend([hnd_xistart, hnd_xiend], 'Location', 'best'); 

% Receiver function. NO ysfhit
axes(ax1); 
set(gca, 'XGrid', 'on', 'XMinorTick', 'on'); 
xlim([-1, 25]); % xlim in seconds. 
yshift_const = 0.00; % Make sure receiver functions aren't shifted up or down based on their index. 

txt_t_xi = '$t(H_{\xi}, \kappa_{\xi}, \xi)$'; % Text for legends in phase timing plot. 
txt_t_0b = '$t(H_{1}, \kappa_{1}, 1)$'; 
txt_t_00 = '$t(H_{\xi}, \kappa_{\xi}, 1)$'; 

nrf = size(rf,2); 
for irf = 1:nrf % Loop over and plot each receiver funcion. 
    yshift = irf * yshift_const; 
    rfi = rf(:,irf);
    hnd_rf = plot(tt, yshift+rfi, 'linewidth', 0.75,...
        'Color', [0,0,0,0.14]); ; 
end

for irf = 1:nrf % Loop over each receiver function and plot the timing where we might predict the Ps etc. pulses. 
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

ylim([-0.25, 0.35]); % Ylimit in amplitude for receiver functions plot. 
lgd = legend([hnd_t_xi, hnd_t_0b, hnd_t_00], 'Location', 'northeast',...
    'Interpreter','latex', 'fontsize', 12); 

exportgraphics(gcf, sprintf('figs/rfs_paper_merged_%s.pdf',sta_name), ...
    'ContentType', 'vector'); % Save figure. 


%% Nice looking HK stack
figure(302); clf; hold on; set(gcf, 'pos', [2623 611 325 212]); 
subplot(1,1,1); hold on; 
set(gca,'ydir', 'reverse', 'LineWidth', 1.5);
box on; 

xlabel('\kappa'); 
ylabel('H (km)'); 
title(sprintf('H-\\kappa stack %s', sta_name), 'fontweight', 'normal'); 

[~,contH] = contourf(gca, K, H, Exi, 15, 'linestyle', 'none'); % Contour of HK energy stack 

cbar = colorbar(gca,'eastoutside'); 
set(cbar, 'fontsize', 12); 

try % Colormap might not be available on some peoples computers. 
    colormap(viridis()); 
catch 
    warning('Missing colormap viridis. Should be in repositories somewhere. '); 
end

scatter(kmax, hmax, 150, 'r', 'pentagram'); 

exportgraphics(gcf, sprintf('figs/hk_nice_%s.jpeg',sta_name),...
    'Resolution', 700); 

%% Receiver function plot, with y shift added based on index. 
% Similar to above receiver function plot, but with vertical shift. 
figure(202); clf; hold on; 
set(gcf, 'pos', [2032 336 578 1029]); 
set(gca, 'LineWidth', 1.5, 'XGrid', 'on', 'XMinorTick', 'on'); box on; 
xlabel('Time (s)'); 
title('Phase timing', 'FontWeight','normal'); 
set(gca, 'YTick', []); 
xlim([-1, 25])
yshift_const = 0.05; % Shift receiver functions by this much times their index. 

nrf = size(rf,2); 
for irf = 1:nrf % Loop over each receiver function and plot it with predicted pulse timings. 
    yshift = irf * yshift_const; 
    rfi = rf(:,irf); % the ith receiver function time series. 

    t_plot_xi = t_pred_xi(irf,:) ; 
    t_plot_00 = t_pred_00(irf,:)'; 
    t_plot_0b = t_pred_0b(irf,:)'; 

    % Scatters are for pulse timing. 
    hnd_t_xi = scatter(...
        t_plot_xi', yshift + interp1(tt, rf(:,irf), t_plot_xi, 'linear'),...
        50, '+', 'CData', c_t_with_xi, 'LineWidth', 2, 'DisplayName', txt_t_xi); % If using true parameters and anisotropic stack
    hnd_t_00 = scatter(...
        t_plot_00', yshift + interp1(tt, rf(:,irf), t_plot_00, 'linear'),...
        50, '+', 'CData', c_t_no_xi, 'LineWidth', 2, 'DisplayName', txt_t_00); % If using true parameters and isotropic stack
    hnd_t_0b = scatter(...
        t_plot_0b', yshift + interp1(tt, rf(:,irf), t_plot_0b, 'linear'),...
        50, '+', 'CData', c_t_iso, 'LineWidth', 2, 'DisplayName', txt_t_0b); % If using parameters from optimizing the isotropic stack

    % The receiver function
    hnd_rf = plot(tt, yshift+rfi, 'k', 'linewidth', 0.75);
end

ylim([-10*yshift_const, yshift_const * (nrf+10)]); % Find a nice looking y limit. Can play with these values. 
lgd = legend([hnd_t_xi, hnd_t_0b, hnd_t_00], 'Location', 'southeast',...
    'Interpreter','latex', 'fontsize', 13); 

text(1, -3*yshift_const     , sprintf('$p=%1.2f$',rayp(1  )), 'interpreter', 'latex' ); % State ray parameter. If it's seconds per degree, should have something like 4.5 to 9. Seconds per kilometer should be something like 0.04 to 0.08 (See Zhu and Kanamori 2000). 
text(1, (nrf+5)*yshift_const, sprintf('$p=%1.2f$',rayp(end)), 'interpreter', 'latex'  ); 

exportgraphics(gcf, sprintf('figs/rfs_%s.pdf',sta_name), ...
    'ContentType', 'vector'); % save figure. 