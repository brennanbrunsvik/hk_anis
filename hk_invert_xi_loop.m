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

% Set some parameters. 
kNum = 51; 
hNum = 50; 
xiNum = 31; 
plot_true = true; 
% % % VP = 6.276; % ELK assumed is 6.276 in IRIS EARS
% % % VS = 1/1.77 * VP; % Specifying because I don't remember which we hold constant. 
% % % zmoh = 30; % ELK estimated at 30 km. Use this for some rough approximations, not so important. 
% % % rhoest = 2.83; 
VP = 6.276; %LB BMN values, mostly the same as for US ELK
VS = 1/1.76 * VP; % Specifying because I don't remember which we hold constant. 
zmoh = 30; % ELK estimated at 30 km. Use this for some rough approximations, not so important. 
rhoest = 2.83; 
acparms = [0.85, 1; 0.5, 30]; 

% Load stations of interest. 
f_sta_list = './models/download_rfs/stationList.csv'; 
sta_name_all = string(table2array(readtable(f_sta_list))); 
sta_name_all = sta_name_all(:,1) + '.' + sta_name_all(:,2); 

sta_name_all = "AE.X16A"; warning('Only doing AE.X16A')

for ista = 1:size(sta_name_all, 1); 
    fprintf('%3.0f sta', ista); 
    sta_name = sta_name_all(ista); 
    
    f = load("models/download_rfs/Ears/gauss_2.5/"+sta_name+"/rfArr.mat"); 
    rf = f.rf; tt = f.tt; rayParmSecDeg = f.rayParmSecDeg; incAngP = f.incAngP; 
    lat = f.stla; lon = f.stlo; 

    if size(rf, 2) < 50; 
        continue; 
    end
    
    rayp = rayParmSecDeg; 
    
    % %% Prepare model. 
    % load('models/model_ta_kmsc.mat'); 
    % sta_name = 'TA.KMSC'; 
    
    % Make a fake simple 2d model. Just for earth flattening. 
    z = 1:100; 
    z = z'; 
    nz = length(z); 
    sz = size(z); 
    rho = rhoest * ones(sz); 
    vs = VS * ones(sz); 
    vp = VP * ones(sz); 
    
    
    % Remove bad receiver functions. 
    % Remove those that have small zero-lag cross-correlation to other receiver functions. 
    
    % auto correlation parameters. These don't need to be perfect. It's just removing really bad receiver functions. 
    % First row: Handle bad parent pulse. Second is for whole receiver function. 
    % First collumn: Cutoff for noramlized average cross correlation to other
    % receiver functions. Second collumn: Time window (absolute value)
    if plot_true; 
        figure(100); clf; hold on; set(gcf, 'pos', [1271 539 933 798]); 
        tiledlayout(2,1,"TileSpacing",'compact'); 
        % % % nexttile(); hold on; box on; xlim([-3, 50]); 
        % % % plot(tt, rf); 
        % % % title(sprintf('%1.0f Receiver functions to start', size(rf,2)), ...
        % % %     'FontWeight','normal'); 
    end
    
    % acparms = [0.93, 1; 0.55, 30]; 
    
    rf_delete = zeros(size(rf,1), 0); 
    for iacremove = 1:2; % Loop over checking the first few seconds of the receiver function, and the whole receiver function. 
        if size(rf,2) < 5; 
            sprintf('Not enough receiver functions keps for %s',sta_name)
            continue
        end
        nrf = size(rf,2); % Number of receiver functions. Changes with each iacremove
        ac = nan(nrf, nrf); % Autocorrelation matrix. 
        twin = abs(tt) < acparms(iacremove,2); % Window within which is the parent pulse, to cross-correlate.  
        ac_cutoff = acparms(iacremove,1); % Get parameters for this iacremove. 
    
        rf_twin = rf(twin,:); 
        rf_twin = rf_twin ./ sqrt(sum(rf_twin .^ 2, 1)); % Normalize
        ac = rf_twin' * rf_twin; % zero-lag cross-correlation
    
        acs = sum(ac) ./ length(ac); % Cross-correlation sum. Averaged. 
        ac_keep = acs > ac_cutoff; % Only keep if there is high-enough average cross-correlation. 
        
        rf_delete = [rf_delete, rf(:,~ac_keep)]; 
        rayp     = rayp(:,ac_keep     ); % Remove ray parameter for receiver functions we aren't keeping. 
        rf       = rf  (:,ac_keep     ); % Remove bad receiver functions. 

    end
    if size(rf,2) < 5; 
        sprintf('Not enough receiver functions keps for %s',sta_name)
        continue
    end
    
    if plot_true; 
        nexttile(); hold on; box on; xlim([-3, 50]); 
        plot(tt, rf, 'blue'); 
        plot(tt, rf_delete, 'red', 'linewidth', 1); 
        title(sprintf('%1.0f Receiver functions kept, %1.0f deleted', size(rf,2), size(rf_delete,2)), ...
            'FontWeight','normal'); 
        
        exportgraphics(gcf, sprintf('figs/%s_rfs.jpeg',sta_name)...
            , 'Resolution', 250); 
    end
    
    %% Sort by ray parameter. Makes plotting more obvious. 
    [rayp, sortp] = sort(rayp); 
    rf = rf(:,sortp); 
    
    %% Loop over xi, doing HK stack 
    
    XI = linspace(0.60, 1.40, xiNum); 
    
    Exi_all = nan(hNum, kNum, xiNum); 
    
    
    for ixi = 1:length(XI); 
        xi = XI(ixi) * ones(sz); 
        phi = 1./xi; 
        eta = 1 * ones(sz); 
    
    
        % Get average crustal values to use with with Zhu and Kanamori math
        [rhoav, vsav, vpav, xiav, phiav, etaav] =  hk_average_crust(...
            rho, vs, vp, xi, phi, eta, z, zmoh); 
    
        xiav = XI(ixi); 
        phiav = 1./xiav; 
        etaav = 1; 
    
        [Exi, H, K, t_pred_all_xi] = hk_anis_loop(...
            rf, tt, rayp, vsav, rhoav, xiav, phiav, etaav, ...
            'kNum', kNum, 'hNum', hNum, 'normalize', false, ...
            'hBounds', [10, 60]); % accounting for xi
    
        Exi_all(:, :, ixi) = Exi; 
    end
    
    e_max = max(Exi_all, [], 'all'); 
    Exi_all = Exi_all / e_max; 
    
    if sum(Exi_all == 1, 'all') > 1; 
        error('Need to get precise indecies of Exi_all max'); 
    else
        [ihmax, ikmax, iximax] = ind2sub(size(Exi_all), find(Exi_all == 1) ); 
        hbest = H(ihmax); 
        kbest = K(ikmax); 
        xibest=XI(iximax); 
        save(sprintf('results/%s.mat',sta_name), 'Exi', 'Exi_all', ...
            'H', 'K', 'XI',...
            'hbest', 'kbest', 'xibest', 'lon', 'lat'); 
    end
    
    %%
    if plot_true; 
        figure(1); clf; hold on; set(gcf, 'pos', [562 747 1089 316]); 
        tiledlayout(1,3,'TileSpacing','compact'); 
        sgtitle(sta_name); 
        cnt_lvl = 0.6:0.05:0.95; 
        
        nexttile(); hold on; box on; 
        exi = Exi_all(ihmax, :, :); 
        exi = reshape(exi, kNum, xiNum); 
        contourf(XI, K, exi, cnt_lvl); 
        plot(XI(iximax), K(ikmax), '^'); 
        xlabel('\xi'); ylabel('\kappa'); 
        
        nexttile(); hold on; box on; 
        exi = Exi_all(:, ikmax, :); 
        exi = reshape(exi, hNum, xiNum); 
        contourf(XI, H, exi, cnt_lvl); 
        plot(XI(iximax), H(ihmax), '^'); 
        xlabel('\xi'); ylabel('H (km)'); 
        set(gca, 'ydir', 'reverse'); 
        
        nexttile(); hold on; box on; 
        exi = Exi_all(:, :, iximax); 
        exi = reshape(exi, hNum, kNum); 
        contourf(K, H, exi, cnt_lvl); 
        plot(K(ikmax), H(ihmax), '^'); 
        xlabel('\kappa'); ylabel('H (km)'); 
        set(gca, 'ydir', 'reverse'); 
%         cmp = flipud(jet(length(cnt_lvl)-1 )); 
%         colormap(cmp); 
%         parula()
        %   
%         cmp = parula(); 
        cbr = colorbar(); 
%         set(cbr.Label.String, 'E/E_{max}')
        cbr.Label.String = 'E/E_{max}'; 
        colormap(viridis()); 

%         caxis([min(cnt_lvl), max(cnt_lvl)]); 
        
        exportgraphics(gcf, sprintf('figs/%s_tradeoff.jpeg',sta_name)...
            , 'Resolution', 200); 
    end
end