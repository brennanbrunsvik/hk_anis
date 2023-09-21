% If you succesfully ran HK stacks on many stations, this is a starting
% script to show some results at many of those stations. You will need to
% change some paths. 

fstainfo ='/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/results_compiled/compiled_results_standard.mat'; % This .mat file has all stations with results used to make 3D mcmc results. Generated from script: /Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/a1_collate_sta_models.m
addpath('/Users/brennanbrunsvik/MATLAB/m_map'); % m_map for plotting. 
addpath('/Users/brennanbrunsvik/MATLAB/borders'); % For getting state borders. 
addpath('./misc')
sta_summary = load(fstainfo).mdls; 
sta_name_all = string(sta_summary.nwk) + "." + string(sta_summary.sta); 

summary = struct(); 

nrf_min = 0; 

% Some QC
hBounds = [20, 50]; % Don't accept if out of this range. 
kBounds = [1.6, 2.1]; 
sta_name_ignore = ["TA.N48A", "TA.L61A"]; 
% sta_name_ignore = ["ASDFASDF"]; 
stacount = 0; 
reason_skip = zeros(length(sta_name_all),1); 
bounds_skip = 0; 
for ista=1:length(sta_name_all); 
    sta_name = sta_name_all(ista); 
    try
        hk_summary = load(sprintf('results/%s/hk_summary.mat', sta_name)); 
    catch
        continue
        fprintf('Could not load results for station %s', sta_name)
    end

    % Manually skip some stations. 
    if any(sta_name == sta_name_ignore)
        disp('Skipping: ignore list')
    end
    skip = any(sta_name == sta_name_ignore); 

    % Correct run version. 
    skip = skip | (~strcmp(hk_summary.run_version, 'version_1')); 
    if (~strcmp(hk_summary.run_version, 'version_1'))
        disp('Skip: wrong version')
    end

    % Have a minimum number of receiver functions. 
    skip = skip | (hk_summary.n_rf < nrf_min); 
    if (hk_summary.n_rf < nrf_min)
        disp('Skip: nrf min')
    end


    % Check for out of bounds
    skip = skip | ...
        any([hk_summary.hmax, hk_summary.hmaxiso] < hBounds(1)) | ...
        any([hk_summary.hmax, hk_summary.hmaxiso] > hBounds(2)) | ...;  
        any([hk_summary.kmax, hk_summary.kmaxiso] < kBounds(1)) | ...;  
        any([hk_summary.kmax, hk_summary.kmaxiso] > kBounds(2));  
%     if any([hk_summary.hmax, hk_summary.hmaxiso] < hBounds(1)) | ...
%         any([hk_summary.hmax, hk_summary.hmaxiso] > hBounds(2)) | ...;  
%         any([hk_summary.kmax, hk_summary.kmaxiso] < kBounds(1)) | ...;  
%         any([hk_summary.kmax, hk_summary.kmaxiso] > kBounds(2))
%         fprintf('Skip due to bounds %s. H: %3.3f, K: %3.3f\n', sta_name, hk_summary.hmax, hk_summary.kmax )
%         bounds_skip = bounds_skip + 1; 
%     end

    % If h or k changed too much, then jumped between local minima. Ignore that result. 
    skip = skip | ...
        (15 < abs((hk_summary.hmax - hk_summary.hmaxiso) / hk_summary.hmax * 100)) | ...
        (15 < abs((hk_summary.kmax - hk_summary.kmaxiso) / hk_summary.kmax * 100)); 
    if (15 < abs((hk_summary.hmax - hk_summary.hmaxiso) / hk_summary.hmax * 100)) | ...
        (15 < abs((hk_summary.kmax - hk_summary.kmaxiso) / hk_summary.kmax * 100))
        disp('Skip: Large change. ') 
    end

    if skip; 
        continue; 
    end

    % We like this result. Continue onward. 
    stacount = stacount + 1; 


    summary(stacount).E00_max     = hk_summary.E00_max; 
    summary(stacount).Exi_max     = hk_summary.Exi_max; 
    summary(stacount).hmax        = hk_summary.hmax   ; 
    summary(stacount).hmaxiso     = hk_summary.hmaxiso; 
    summary(stacount).kmax        = hk_summary.kmax   ; 
    summary(stacount).kmaxiso     = hk_summary.kmaxiso; 
    summary(stacount).n_rf        = hk_summary.n_rf   ; 
    summary(stacount).xi          = hk_summary.xiav   ; 
    summary(stacount).run_version = hk_summary.run_version; 
    summary(stacount).lat         = sta_summary.lat(ista) ; 
    summary(stacount).lon         = sta_summary.lon(ista) ; 
    summary(stacount).sta_name    = sta_name              ;

end
disp('Number skipped due to being out of bounds: ')
disp(bounds_skip)

%%

scatter_size = 25; 

lon = [summary.lon]'; 
lat = [summary.lat]'; 
h   = [summary.hmax]'; 
hiso = [summary.hmaxiso]'; 
xi = [summary.xi]'; 
sta_name = [summary.sta_name]'; 

% keep = abs(xi-1) > .02; 
% lon = lon(keep); 
% lat = lat(keep); 
% h = h(keep); 
% hiso = hiso(keep);  
% xi = xi(keep); 

dh_perc = (h-hiso)./h*100; 

% lon = lon(xi>1.05); 
% lat = lat(xi>1.05); 
% h = h(xi>1.05); 
% hiso = hiso(xi>1.05); 
% xi = xi(xi>1.05); 


% Plots
clims = [20, 60]; 
lolim = [-87, -71]; 
lalim = [27 ,  46]; 
figure(1); clf; hold on; set(gcf, 'pos', [297 391 495 486]); 
tiledlayout(2, 2, 'TileSpacing','compact'); 

letter_shift = 1.07; 
letter_shift_left = -0.14;
nexttile(); hold on; 
text(letter_shift_left, letter_shift, '(a)', 'units', 'normalized', 'fontsize', 12)
[stax, stay] = fun_mapbackground(lolim, lalim, lon, lat); 
scatter(stax, stay , scatter_size, h, 'filled'); 
colorbar(); 
clim(clims); 
title('H (km)', 'fontweight', 'normal');

nexttile(); hold on; 
text(letter_shift_left, letter_shift, '(b)', 'units', 'normalized', 'fontsize', 12)
[stax, stay] = fun_mapbackground(lolim, lalim, lon, lat); 
scatter(stax, stay, scatter_size, xi, 'filled');
% scatter(lon, lat , scatter_size, hiso, 'filled');
colorbar(); 
% clim(clims); 
title('\xi', 'fontweight', 'normal');

nexttile(); hold on; 
text(letter_shift_left, letter_shift, '(c)', 'units', 'normalized', 'fontsize', 12)
[stax, stay] = fun_mapbackground(lolim, lalim, lon, lat); 
scatter(stax, stay, scatter_size, dh_perc, 'filled');
colorbar(); 
title('Change in H, percent', 'fontweight', 'normal');


nexttile(); hold on; 
text(letter_shift_left, letter_shift, '(d)', 'units', 'normalized', 'fontsize', 12)
title('Percent change in H versus \xi', 'fontweight', 'normal');
scatter(xi, dh_perc); 
% text(xi, dh_perc, sta_name); % Use these to see which stations are weird 
xlabel('\xi'); 
ylabel('Change in H, percent')
box on; grid on; set(gca, 'linewidth', 1.5);

dhp_xi = dh_perc'/(xi-1)'; % The slope of dh to (xi-1). Do a simple least squares matrix fit, where some number * (xi-1) gives h. Use matlab \ to do least squares solver. 
dh_pred = dhp_xi .* (xi-1); 
plot(xi, dh_pred, 'linewidth', 1.5); 

fprintf('Solved slope: change from 1 to 1.01 and change the h solution percent by %3.3f\n', dhp_xi/100)

set(gcf, 'Renderer', 'painters'); 
exportgraphics(gcf, 'figs/hk_from_mcmc_models.pdf', 'ContentType', 'vector'); 
