% IN PROGRESS, not yet documented. 
% Make plots to show analytical error from ignoring anisotropy. 

clear; 
clc; 
restoredefaultpath; 

addpath("elastic/"); 
addpath("hk_calculation/"); 
addpath('~/Documents/repositories/Base_code/colormaps/colormaps_ander_biguri'); 

%%rayp
% rayp = linspace(4, 9, 100); 
% rayp = 8; 
raypNum = 17; 
kNum = 15; 
hNum = 16; 

rayp_all = linspace(4, 9, raypNum); 
rayp_conv = @(r)r/(6371 * 2 * pi/360); 
vs = 3.8 .* ones(1);
% vp = 1.75 .* ones(size(rayp));
rho = 2.6 .* ones(1);
xi = 1.15; 



p_change         = zeros(kNum, raypNum); 
p_change_t_ps    = zeros(kNum, raypNum); 
p_change_t_ppss  = zeros(kNum, raypNum); 
p_change_t_ppps  = zeros(kNum, raypNum); 

%%
for irayp = 1:raypNum; 
    rayp = rayp_all(irayp); 
    [p_change(:,irayp), p_change_t_ps(:,irayp), ...
        p_change_t_ppss(:,irayp), p_change_t_ppps(:,irayp), kVec, hVec] = ...
        hk_anis_error_calc(rayp, vs, rho, xi, 1./xi, 1, 'kNum', kNum, 'hNum', hNum); 
end

%%


clim_man = [0, 7]; 
% clim = 
% reds = [0:max(clim_man)]./(max(clim_man)); 
% cmap_cust = [reds; zeros(size(reds)); zeros(size(reds)) ]'; 
% reds = [1:max(clim_man)]./(max(clim_man)); 
reds = ones(1, max(clim_man));  
greens = linspace(1,0, length(reds)); 
blues = linspace(1,0,length(reds)); 
cmap_cust = [reds; greens; blues ]'; 
cmap_cust = jet(max(clim_man) .* 2 ); 
clines = [0:max(clim_man)]; 

toabs = @(f)abs(f); 

figure(11); clf; hold on; 
tiledlayout(2,2, 'TileSpacing','compact'); 
sgtitle('Timing error'); 

nexttile(); box on; set(gca, 'LineWidth', 1.5); caxis(clim_man); 
contourf(kVec', rayp_all, toabs(p_change'), clines); 
caxis(clim_man); 
colormap(cmap_cust); 
ylabel('p (s/km)'); 
title('Average', 'fontweight', 'normal'); 

nexttile(); box on; set(gca, 'LineWidth', 1.5);  caxis(clim_man); 
contourf(kVec', rayp_all, toabs(p_change_t_ps'), clines); 
caxis(clim_man); 
title('Ps', 'fontweight', 'normal');  


nexttile(); box on; set(gca, 'LineWidth', 1.5);  caxis(clim_man); 
contourf(kVec', rayp_all, toabs(p_change_t_ppss'), clines);  
caxis(clim_man); 
xlabel('\kappa'); 
ylabel('p (s/km)'); 
title('PpSs, PsPs', 'fontweight', 'normal'); 


nexttile(); box on; set(gca, 'LineWidth', 1.5);  
contourf(kVec', rayp_all, toabs(p_change_t_ppps'), clines);  
caxis(clim_man); 
xlabel('\kappa'); 
title('PpPs', 'fontweight', 'normal'); 

cbar = colorbar(); 
set(cbar.Label, 'String', 'H error percent'); 
cbar_pos = cbar.Position; 
cbar_pos(4) = 0.5; 
% cbar_pos(2) = 0.33; 
cbar_pos(2) = 0.25; 
cbar_pos(1) = cbar_pos(1) + 0.01; 
set(cbar, 'Position', cbar_pos); 

exportgraphics(gcf, sprintf('figs/timing_error_xi%1.2f.pdf', xi)); 


