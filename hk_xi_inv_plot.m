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
addpath('/Users/brennanbrunsvik/MATLAB/m_map'); 
addpath('/Users/brennanbrunsvik/MATLAB/borders'); 
addpath('/Users/brennanbrunsvik/Documents/repositories/Base_code/colormaps/redblue'); 

path_results = 'results/'; 

files = string(split(ls(path_results)));
files = files(1:end-1); % Empty last value, for some reason. 

%%

lon = []; 
lat = []; 
xi = []; 
anis_improvement = []; 

for ifile = 1:length(files); 
    file = files(ifile);
    f = load("results/"+file); 
    try 
        H = f.H; 
        K = f.K; 
        XI = f.XI; 
        Exi = f.Exi; 
        hbest = f.hbest; 
        kbest = f.kbest; 
        xibest = f.xibest; 
        Exi_all = f.Exi_all; 
    catch
        continue
    end

    nh = size(Exi_all,1); 
    nk = size(Exi_all,2); 
    nxi = size(Exi_all,3); 

%     Exi_iso = Exi_all(); 

    % Temporary, get best E from iso and aniso stacks
    Exi_p = permute(Exi_all, [3,1,2]); 
    Exi_1 = interp1(XI', Exi_p, 1, 'spline'); 
    Exi_1 = reshape(Exi_1, size(Exi_1,2), size(Exi_1,3)); 

    Exi_p = permute(Exi_all, [3,1,2]); 
    Exi_b = interp1(XI', Exi_p, xibest, 'spline'); 
    Exi_b = reshape(Exi_b, size(Exi_b,2), size(Exi_b,3)); 

    Ebest_1  = max(Exi_1, [], 'all'); 
    Ebest_xi = max(Exi_b, [], 'all'); 

    anis_improvement_i = Ebest_xi ./ Ebest_1; 

    if anis_improvement_i < 1.05; 
        continue; 
    end

    lon = [lon, f.lon]; 
    lat = [lat, f.lat]; 
    xi  = [xi , f.xibest]; 
    anis_improvement = [anis_improvement, anis_improvement_i]; 


%     contour4(Exi_all, 0.9); 

end


cvals = [[0.6:0.05:0.95]];% .* max(max(Exi)); 
figure(103); clf; hold on; 
tiledlayout(1,2,'TileSpacing','compact'); 
nexttile(); hold on; 
contourf(Exi_1, cvals); 
nexttile(); hold on; 
contourf(Exi_b, cvals); 




figure(1); clf; hold on; box on; set(gca, 'linewidth', 1.5); 
m_proj('miller','long',[min(lon)-1 max(lon)+1],'lat',[min(lat)-1, max(lat+1)]);
m_coast('patch',[1 1 1]); % m_coast('patch',[1 .85 .7]);
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
end
m_coast('patch',[1 .85 .7]);
% m_elev('contourf',[500:500:6000]);
m_coast(); 
m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',[.3 .75 1]);
% m_elev('contour');


[stax, stay] = m_ll2xy(lon, lat);  
scatter(stax, stay, 60, anis_improvement, 'filled');  


% colormap('cool'); 
colormap(viridis()); 
cbar = colorbar(); cbar.Label.String = 'E_{\xi}/E_{\xi=1}'; 

% m_shadedrelief(); 
% m_elev('contourf',[500:500:6000], 'image');



title('\xi Inverted', 'fontweight', 'normal'); 
exportgraphics(gcf, 'figs/xi_inv_results.pdf'); 


%%
    

figure(2); clf; hold on; box on; set(gca, 'linewidth', 1.5); 
m_proj('miller','long',[min(lon)-1 max(lon)+1],'lat',[min(lat)-1, max(lat+1)]);
m_coast('patch',[1 1 1]); % m_coast('patch',[1 .85 .7]);
[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
end
m_coast('patch',[1 .85 .7]);
% m_elev('contourf',[500:500:6000]);
m_coast(); 
m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',[.3 .75 1]);
% m_elev('contour');


[stax, stay] = m_ll2xy(lon, lat);  
scatter(stax, stay, 60, xi, 'filled');  


% colormap('cool'); 
colormap('redblue'); 
caxis([min(XI), max(XI)]); 
cbar = colorbar(); cbar.Label.String = '\xi'; 

% m_shadedrelief(); 
% m_elev('contourf',[500:500:6000], 'image');



title('\xi Inverted', 'fontweight', 'normal'); 
exportgraphics(gcf, 'figs/xi_inv_results.pdf'); 
% f = load("models/download_rfs/Ears/gauss_2.5/"+sta_name+"/rfArr.mat"); 
% rf = f.rf; tt = f.tt; rayParmSecDeg = f.rayParmSecDeg; incAngP = f.incAngP; 
% lat = f.stla; lon = f.stlo; 

%% Second try at figure
figure(103); clf; hold on; set(gcf, 'pos', [1057 948 440 543]); 
[A,R,ATTRIB] = readBasemapImage("landcover", [min(lat), max(lat)], [min(lon), max(lon)]); 
[x,y] = projfwd(R.ProjectedCRS,lat,lon);
mapshow(A,R); 
axis off; 
scatter(x, y, 60, xi, 'filled'); 
colormap('redblue'); 
caxis([min(XI), max(XI)]); 
cbar = colorbar(); cbar.Label.String = '\xi'; 
exportgraphics(gcf, 'figs/xi_inv_results_v2.pdf'); 
