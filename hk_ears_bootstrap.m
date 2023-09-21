% This script does bootstrap resampling of HK stacks. It will need some
% modification in terms of providing the velocity and xi models, and
% changing paths, to run on your computer. 

clear; 
clc; 
restoredefaultpath; 

addpath("elastic/"); 
addpath("hk_calculation/"); 
addpath("hk_calculation/"); 
addpath('~/Documents/repositories/Base_code/colormaps/colormaps_ander_biguri'); 

% fstainfo contains list of stations used for MCMC 2023. Can make a different station list if you want. 
fstainfo ='/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/results_compiled/compiled_results_standard.mat'; % This .mat file has all stations with results used to make 3D mcmc results. Generated from script: /Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/a1_collate_sta_models.m

vs_model_style = "mcmc2023"; % mcmc2023 to get absolute vs, xi, etc. from Brunsvik et al. submitted (2023?). MCMC paper. Alternatively, Define vs, etc., manually, and get receiver functions from ears (see hk_ears_onesta.m for example of how this works). 

run_version = 'version_1'; % Starting versions 2023/07/08

kNum = 200; 
hNum = 201; 
plot_true = true; 

% Limit the bounds for this script, because the bootstrapping is slow. 
hBounds = [20, 55]; 
kBounds = [1.6, 2.1]; warning('changing bounds')

nrfmin = 30; 
acparms = [0.85, 1; 0.45, 40]; % Might be more gentle here when working with many receiver functions, or most receiver functions will not be used. 


eta_single = 1; 
phi_xi_const = 1; % define 1/phi = phi_xi_const * xi. Still should double check how this is implimented in code below, brb2023.03.21

load_from_orig = true; % Set false for TA.KMSC if you don't have access to the original MCMC outputs. Load from original data formats. ~100 Mb per station. Can condense the model outputs and save only the key parameters, but still will take about 3-4 Gb for 400 stations. I include a data file for one station so that a user can try on that. 
save_condensed_data = false; % Save a file with just needed parameters for HK. ~9 Mb/station. 
rm_old_fold = true; if rm_old_fold; warning('Removing previous station results. '); end 
show_all_plots = false; 
% hi_resolution_figs = false; 


%% Prepare model. 
% sta_name = 'AE.X16A'; 

if ~ (vs_model_style == "mcmc2023"); 
    error('TODO Uncomment below, and debug')
    % % Load stations of interest. 
    % f_sta_list = './models/download_rfs/stationList.csv'; 
    % sta_name_all = string(table2array(readtable(f_sta_list))); 
    % sta_name_all = sta_name_all(:,1) + '.' + sta_name_all(:,2); 
elseif vs_model_style == "mcmc2023"; 
%     sta_name_all = ["TA.KMSC"]; 
    slist = load(fstainfo).mdls; 
    sta_name_all = string(slist.nwk) + "." + string(slist.sta); 
end



sta_name_all = ["TA.KMSC"]'; % If you only want to run a few stations for new, replace sta name with those here. 
hi_res_stas = ["US.GOGA", "TA.KMSC"]'; % Make higher resolution figures for these stations only. 
for ista = 1:length(sta_name_all); 

sta_name = sta_name_all(ista); 

res_export = 200; % Choose base figure resolution
if any(sta_name==hi_res_stas); % Some stations can have higher figure resolution
    res_export = res_export + 500; 
end

fprintf("\n ----- %s ---- nsta = %5.0f / %5.0f \n", sta_name, ista, length(sta_name_all))

if rm_old_fold & (exist(sprintf('./results/%s', sta_name))~=0); 
    rmdir(sprintf('./results/%s', sta_name), 's'); 
end

% Load model. Can use mcmc format results. Or get results from ears. Both require some manual labor to get model files ready. 
if ~ (vs_model_style == "mcmc2023"); % Not using our mcmc results? define your own vs, xi, etc. values, and get receiver functions from ears. 
    try 
        % Put your appropriate path here. 
        f = load("/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/repositories/hk_anis/models/download_rfs/Ears/gauss_2.5/"+sta_name+"/rfArr.mat"); 
    catch
        fprintf('Tried loading from extDrive. Did not work. Trying from models/download_rfs/Ears/...\n')
        f = load("models/download_rfs/Ears/gauss_2.5/"+sta_name+"/rfArr.mat"); 
    end
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

    z = 1:100; 
    z = z'; 
    xi = ones(size(z)) .* xi; 
    eta = ones(size(z)) .* eta_single; 
    phi = 1./phi_xi_const .* 1./xi; % See paper for description of this relationship. 
    nz = length(z); 
    sz = size(z); 
    rho = rhoest * ones(sz); 
    vs = VS * ones(sz); 
    vp = VP * ones(sz); 
    
    % Get average crustal values to use with with Zhu and Kanamori math
    [rhoav, vsav, vpav, xiav, phiav, etaav] =  hk_average_crust(...
        rho, vs, vp, xi, phi, eta, z, zmoh); 
else % Load mcmc style model. 
    % Should load something like this directly into matlab: 
        % %      eta: [245×1 double]
        % %      phi: [245×1 double]
        % %     rayp: [6.4886 5.1054 7.3144 4.9234 8.7014 5.6641 6.1073 8.2694 … ]
        % %       rf: [8190×164 double]
        % %      rho: [245×1 double]
        % %       tt: [8190×1 double]
        % %       vp: [245×1 double]
        % %     vpvs: 1.9346
        % %       vs: [245×1 double]
        % %       xi: [245×1 double]
        % %        z: [245×1 double]
        % %     zmoh: 29.6104

    sta_name_split = sta_name.split('.');
    sta = sta_name_split(2); 
    net = sta_name_split(1); 

    if load_from_orig; 
        % TODO  a functino to automatically get mcmc result and make the .mat file. 
        model_fold_orig = "/Volumes/extDrive/offload/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/STASinv_collate"% Determine path. 
        stamp = 'standard'; 
    
        modin = load(sprintf("%s/%s_%s_dat1/%s/final_model.mat", ...
            model_fold_orig, sta, net, stamp)).final_model; 
        datin = load(sprintf("%s/%s_%s_dat1/%s/final_predata.mat", ...
            model_fold_orig, sta, net, stamp)).final_predata; 
        fmodout = sprintf("./models/model_%s_%s.mat", lower(net), lower(sta) ); % Where to save these results to. 
    
        % Load from model. 
        zmoh = modin.zmohav; 
        vpvs = modin.vpvsav; 
        xi = 1+modin.Sanisav/100; % brb2023/07/08. Verified this gives same results as what I had before first paper submission for TA.KMSC. For the MCMC paper, percent anisotropy was simply xi*100 - 100, for convenience/code compatibility. 
        z = modin.Z; 
        vs = modin.VSav; 
        vp = modin.VPav; 
        rho = modin.rhoav; 
    
        % Derive phi and eta how we want for this paper. Could optionally take phi and eta from MCMC inverted model. 
        phi = 1./phi_xi_const .* 1./xi; % Not the loaded models phi, but the phi we want for this paper. 
        eta = ones(size(z)) .* eta_single;  
    
        % Load from data. 
        waves = datin.HKstack_P.waves; 
        rf = waves.rf; 
        tt = waves.tt; 
        rayp = waves.rayParmSecDeg; % brb2023/07/08 I verified the this is the units as used before the papers first submission. Should be fine. Values at about 7. 
    
        % Save model for simplicity in the future. 
        if save_condensed_data; 
            save(fmodout, 'eta', 'phi', 'rayp', 'rf', 'rho', 'tt', 'vp', 'vpvs', 'vs', 'xi', 'z', 'zmoh'); 
        end
    else % Load condensed data file. 
        load(sprintf('models/model_%s_%s.mat', lower(net), lower(sta)) ); % Loads rf, tt, eta, phi, rayp, etc. 
    end

    % Get average crustal values to use with with Zhu and Kanamori math
    [rhoav, vsav, vpav, xiav, phiav, etaav] =  hk_average_crust(...
        rho, vs, vp, xi, phi, eta, z, zmoh); 

    % TODO brunsvik et al submitted (2023) has phi = 1. Could use phi_xi_const to change phi here. 
%     phi_xi_const = 1; % define 1/phi = phi_xi_const * xi. Still should double check how this is implimented in code below, brb2023.03.21

end
 

%% Remove bad receiver functions. 
% Remove those that have small zero-lag cross-correlation to other receiver functions. 

% cross correlation parameters. These don't need to be perfect. It's just removing really bad receiver functions. 
% First row: Handle bad parent pulse. Second is for whole receiver function. 
% First collumn: Cutoff for noramlized average cross correlation to other
% receiver functions. Second collumn: Time window (absolute value)

fprintf("\nCleaning receiver functions with correlation of\n" + ...
    "%1.2f from -%1.2f to -%1.2f seconds and \n%1.2f from -%1.2f to %1.2f seconds\n", ...
    acparms(1,1), acparms(1,2), acparms(1,2), acparms(2,1), acparms(2,2), acparms(2,2))

fprintf('\nStarting with %1.0f receiver functions.\n', size(rf,2))

figure(61); clf; hold on; 
subplot(2,1,1); hold on; 
title('All')
plot(tt, rf); 
xlim([-5, 40]); 

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

n_rf = size(rf,2)

fprintf('\nEnding with %1.0f receiver functions.\n', n_rf)


if n_rf < nrfmin; % Not saving these results. 
    fprintf('Not enough receiver functions for %s (%5.0f/%5.0f). Continue', sta_name, n_rf, nrfmin); 
    continue
else % We are saving these results. Prepare to save them. 
    if ~ exist(sprintf('./results/%s', sta_name)); 
        mkdir(sprintf('./results/%s', sta_name)); 
    end
end

subplot(2,1,2); hold on; 
title('Kept'); 
plot(tt, rf); 
xlim([-5, 40]); 

%% Sort by ray parameter. Makes plotting more obvious. 
[rayp, sortp] = sort(rayp); 
rf = rf(:,sortp); 
rf_orig = rf; % Save original here. We can scale rf in individual cells, if we copy rf from rf_orig. 

%% Loop over receiver functions making HK stack for each. 
% Note: Exi, include xi. E00, standard process without xi.  
nrf = size(rf,2); 
n_bootstraps = 200; 
indx_boots = balanced_resampling(nrf, n_bootstraps); 

h_boots = nan(n_bootstraps,1); 
k_boots = nan(n_bootstraps,1); 


% figure(1); clf; hold on; 
% tiledlayout(5, 5); 
% hBounds = [20, 55]; 
% kBounds = [1.6, 2.1]; warning('changing bounds')
recalc = false; 
if recalc
    parfor iboot = 1:n_bootstraps;  

        rf_resamp = rf(:, indx_boots(:,iboot)); 
        [Exi, H, K, t_pred_all_xi] = hk_anis_loop(...
            rf_resamp, tt, rayp, vsav, rhoav, xiav, phiav, etaav, ...
                    'kNum', kNum, 'hNum', hNum, 'normalize', false, ...
                    'hBounds', hBounds); % accounting for xi
        [ihmax, ikmax] = find(Exi == max(Exi ,[], 'all')); % If you were to optimize while accounting for anisotropy, which k and h have the highest energy?
        kmax = K(ikmax); 
        hmax = H(ihmax); 
        
        h_boots(iboot) = hmax; 
        k_boots(iboot) = kmax; 
        
        % nexttile(); hold on; 
        % contourf(H, K, Exi'); 
        % scatter(hmax, kmax, 100, 'filled'); 
        
        % % Uncomment to show errors. 
        % fprintf('iboot = %1.0f\n',iboot)
        % fprintf('H=%3.6f. K=%3.6f\n', hmax, kmax)
        
    
    end
    save(sprintf('results/boot_res/boot_res_%s.mat', sta_name), 'h_boots', 'k_boots', 'hBounds', 'kBounds' ); 
else 
    load(sprintf('results/boot_res/boot_res_%s.mat', sta_name))
end

%%


%%
fprintf('Figure 1 is x and y bounds set for a very tight standard deviation right now: change if your plots dont make sense\n')

n_hist = 7; 

h_std = std(h_boots); 
k_std = std(k_boots); 
h_mean = std(h_boots); 
k_mean = std(k_boots); 

hplt_bounds = [24, 36]; 
kplt_bounds = [1.8, 2.1]; 

figure(1); clf; hold on; set(gcf, 'pos', [632 116 467 434]); 
tiledlayout(2,2, 'TileSpacing','compact'); 
sgtitle(sprintf('Station %s', sta_name)); 
nexttile(); cla; hold on; set(gca, 'LineWidth', 1.5); box on; grid on; 
histogram(h_boots, n_hist, 'Orientation','horizontal'); 
ylabel('H'); 
xlabel('Count')
ylim(hplt_bounds); 

nexttile(); cla; hold on; set(gca, 'LineWidth', 1.5); box on; grid on; 
scatter(k_boots, h_boots, 10, 'filled'); 
ylabel('H'); 
xlabel('\kappa')
ylim(hplt_bounds); 
xlim(kplt_bounds); 
text(.1, .88, sprintf('std(H)=%1.4f', h_std), 'units', 'normalized'); 
text(.1, .81, sprintf('std(K)=%1.4f', k_std), 'units', 'normalized'); 


nexttile(4); cla; hold on; set(gca, 'LineWidth', 1.5); box on; grid on; 
histogram(k_boots, n_hist); 
xlabel('\kappa'); 
ylabel('Count'); 
xlim(kplt_bounds); 

exportgraphics(gcf, sprintf('results/boot_res/boot_res_%s.pdf', sta_name))

end