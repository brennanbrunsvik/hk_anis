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

% fstainfo contains list of stations used for MCMC 2023. Can make a different station list if you want. 
fstainfo ='/Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/data/results_compiled/compiled_results_standard.mat'; % This .mat file has all stations with results used to make 3D mcmc results. Generated from script: /Users/brennanbrunsvik/Documents/UCSB/ENAM/THBI_ENAM/functions/plotting/many_stas/a1_collate_sta_models.m

vs_model_style = "mcmc2023"; % mcmc2023 to get absolute vs, xi, etc. from Brunsvik et al. submitted (2023?). MCMC paper. Alternatively, Define vs, etc., manually, and get receiver functions from ears (see hk_ears_onesta.m for example of how this works). 

run_version = 'version_1'; % Starting versions 2023/07/08

kNum = 510; 
hNum = 500; 
% % kNum = 200; 
% % hNum = 201; 
% xiNum = 31; 
plot_true = true; 
hBounds = [10, 55];

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



% nrf_used = zeros(size(sta_name_all)); 

tic; 
% sta_name_all = ["US.GOGA", "TA.KMSC"]'; % If you only want to run a few stations for new, replace sta name with those here. 
sta_name_all = ["TA.KMSC"]'; % If you only want to run a few stations for new, replace sta name with those here. 

hi_res_stas = ["US.GOGA", "TA.KMSC"]'; % Make higher resolution figures for these stations only. 
for ista = 1:length(sta_name_all); 

sta_name = sta_name_all(ista); 

res_export = 200; % Choose base figure resolution
if any(sta_name==hi_res_stas); % Some stations can have higher figure resolution
    res_export = res_export + 500; 
end

fprintf("\n ----- %s ---- nsta = %5.0f / %5.0f \n", sta_name, ista, length(sta_name_all))
toc()

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

fprintf('\nHmax_xi = %1.3f, Hmax_iso = %1.3f. dH = %1.3f. percent dH = %1.3f\n', ...
    hmax, hmaxiso, hmax-hmaxiso, (hmax-hmaxiso)/hmaxiso*100)
fprintf('\nKmax_xi = %1.3f, Kmax_iso = %1.3f. dK = %1.3f. percent dK = %1.3f\n', ...
    kmax, kmaxiso, kmax-kmaxiso, (kmax-kmaxiso)/kmaxiso*100)
fprintf('\nRatio of E_xi over E_iso = %1.2f\n', Exi_max/E00_max)


fprintf('\n --- Crustal thinning exercise --- \n\n')
ce_hstart = 40; % Height starting
ce_wstart = 300; % Width starting
ce_wend   = ce_wstart * hmax / ce_hstart; 
ce_wend_i = ce_wstart * hmaxiso / ce_hstart;
fprintf('Width start: %1.3f\n', ce_wstart) 
fprintf('Width end anisotropic: %1.3f\n', ce_wend)
fprintf('Width end isotropic: %1.3f\n', ce_wend_i)
fprintf('Width difference: %1.3f\n', ce_wend-ce_wend_i)
fprintf('Beta anisotropic: %1.3f\n', ce_hstart/hmax)
fprintf('Beta isotropic: %1.3f\n', ce_hstart/hmaxiso)
fprintf('\n --- END Crustal thinning exercise --- \n')


save(sprintf('results/%s/hk_summary.mat', sta_name), 'hmax', 'hmaxiso', 'kmax', 'kmaxiso', 'Exi_max', 'E00_max', 'n_rf', 'run_version', 'vsav', 'rhoav', 'xiav', 'phiav', 'etaav'); 


%% Start figures. 
rf = rf_orig * 2.7; % DELETE THE _orig PART, should not cause a problem

% Some custom colors. 
c_xilow = [179, 54, 214]./255; % Low values of xi
c_xihi = [201, 135, 2]./255; 
c_mod_true = 'y'; % Yellow, star? 
c_t_with_xi = [3, 15, 252]./255; ;% Anisotropic 
c_t_no_xi = [230, 2, 14]./255; 
c_t_iso = [10, 247, 30]./255; 
cnt_mod = 1/1.5; % Multiply colors by this when doing contours. Make them a bit darker. 
xlbl = 0.0; % x (and y) positions for (a) and (b) labelling, for subfigures. 
ylbl = 1.05; 


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
% yticklabels([]); 
ylabel('Ray parameter (km/s)'); 
title('Receiver functions', 'FontWeight','normal', 'Interpreter', 'latex');  
text(xlbl, ylbl, '(a)', 'Units','normalized'); 


axes(ax2); 
xlabel('\kappa'); 
ylabel('H (km)'); % yticklabels([]); 
grid on;
title('$H\kappa$ stack', 'fontweight', 'normal', 'Interpreter', 'latex'); 
text(xlbl, ylbl, '(b)', 'Units','normalized'); 


% HK stack 
axes(ax2); 
fhand_norm = @(inval)inval ./ max(max(inval)); % Return normalized inval 
set(gca,'ydir', 'reverse');

% plt_ylim = [29, 39]; 
% plt_xlim = [1.6, 2.1]; 
plt_ylim = hmax + 5.*[-1, 1]; 
plt_xlim = kmax + .25.*[-1, 1]; 
% Shift limits to within range of H and K. 
% TODO do this for all 4 H, K extremes. 
if max(plt_xlim) > max(K); 
    plt_xlim = plt_xlim - (max(plt_xlim) - max(K)); % Make sure limits aren't out of bounds
end
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
set(gca, 'YGrid', 'on', 'YMinorTick', 'off'); 
xlim([-1, 30]); % xlim in seconds. 
yshift_const = 1; % Make sure receiver functions aren't shifted up or down based on their index. 

txt_t_xi = '$t(H_{\xi}, \kappa_{\xi}, \xi)$'; % Text for legends in phase timing plot. 
txt_t_0b = '$t(H_{1}, \kappa_{1}, 1)$'; 
txt_t_00 = '$t(H_{\xi}, \kappa_{\xi}, 1)$'; 

nrf = size(rf,2); 
yshift_all = rayp' .* yshift_const; 
yshift_all = round(yshift_all); % Bin the receiver functions 


%%% brb2023/07/09 Start plot replacement code. Replaced the slow loop with way faster plotting. TODO
%%% bring this over to the single station plotting script. 
% % % for irf = 1:nrf % Loop over each receiver function and plot the timing where we might predict the Ps etc. pulses. 
% % %     yshift = irf * yshift_const; 
% % %     rfi = rf(:,irf);
% % % 
% % %     t_plot_xi = t_pred_xi(irf,:); 
% % %     t_plot_00 = t_pred_00(irf,:)'; 
% % %     t_plot_0b = t_pred_0b(irf,:)'; 
% % % 
% % %     LWS = 1; 
% % %     sct_symb = '.'; 
% % %     hnd_t_xi = scatter(...
% % %         t_plot_xi', yshift + interp1(tt, rf(:,irf), t_plot_xi, 'linear'),...
% % %         50, sct_symb, 'LineWidth', LWS, 'DisplayName', txt_t_xi,...
% % %         'CData', c_t_with_xi); % If using true parameters and anisotropic stack
% % %     hnd_t_00 = scatter(...
% % %         t_plot_00', yshift + interp1(tt, rf(:,irf), t_plot_00, 'linear'),...
% % %         50, sct_symb, 'LineWidth', LWS, 'DisplayName', txt_t_00,...
% % %         'CData', c_t_no_xi); % If using true parameters and isotropic stack
% % %     hnd_t_0b = scatter(...
% % %         t_plot_0b', yshift + interp1(tt, rf(:,irf), t_plot_0b, 'linear'),...
% % %         50, sct_symb, 'LineWidth', LWS, 'DisplayName', txt_t_0b,...
% % %         'CData', c_t_iso); % If using parameters from optimizing the isotropic stack
% % % end
t_plot_xi_all = zeros(nrf, 3); 
t_plot_00_all = zeros(nrf, 3); 
t_plot_0b_all = zeros(nrf, 3); 
y_plt_xi_all  = zeros(nrf, 3); 
y_plt_00_all  = zeros(nrf, 3); 
y_plt_0b_all  = zeros(nrf, 3); 
for irf = 1:nrf % Loop over each receiver function and plot the timing where we might predict the Ps etc. pulses. 
    yshift = yshift_all(irf); 
    rfi = rf(:,irf);

    t_plot_xi = t_pred_xi(irf,:); 
    t_plot_00 = t_pred_00(irf,:)'; 
    t_plot_0b = t_pred_0b(irf,:)'; 

    t_plot_xi_all(irf,:) = t_plot_xi; 
    t_plot_00_all(irf,:) = t_plot_00; 
    t_plot_0b_all(irf,:) = t_plot_0b; 
    y_plt_xi_all(irf,:)  = yshift + interp1(tt, rf(:,irf), t_plot_xi, 'linear'); 
    y_plt_00_all(irf,:)  = yshift + interp1(tt, rf(:,irf), t_plot_00, 'linear'); 
    y_plt_0b_all(irf,:)  = yshift + interp1(tt, rf(:,irf), t_plot_0b, 'linear'); 

end

% %%% Plot one summed receiver function. 
% rf_sum = sum(rf, 2); 
% rf_sum = rf_sum / max(rf_sum) * max(rf, [], 'all'); % Normalize
% rms_rf_sum = rms(rf_sum); % clip
% clip_rf = 2 * rms_rf_sum; % clip
% rf_sum(rf_sum> clip_rf) =   clip_rf; % clip
% rf_sum(rf_sum<-clip_rf) = - clip_rf; % clip
% yshift_sum = min(yshift_all)-.05*(max(yshift_all)-min(yshift_all)); 
% hnd_rf_sum = plot(tt, yshift_sum+rf_sum, 'linewidth', 2*0.75,...
%         'Color', [0,0,0,1]); 

ylim([min(yshift_all) - .5, max(yshift_all) + 1 ]); 

LWS = 1; 
sct_symb = '.'; 
hnd_t_xi = scatter(...
    t_plot_xi_all, y_plt_xi_all,...
    50, sct_symb, 'LineWidth', LWS, 'DisplayName', txt_t_xi,...
    'CData', c_t_with_xi); % If using true parameters and anisotropic stack
hnd_t_00 = scatter(...
    t_plot_00_all, y_plt_00_all,...
    50, sct_symb, 'LineWidth', LWS, 'DisplayName', txt_t_00,...
    'CData', c_t_no_xi); % If using true parameters and isotropic stack
hnd_t_0b = scatter(...
    t_plot_0b_all, y_plt_0b_all,...
    50, sct_symb, 'LineWidth', LWS, 'DisplayName', txt_t_0b,...
    'CData', c_t_iso); % If using parameters from optimizing the isotropic stack
%%% End plot code replacement 

% rf = rf * 3; % scale receiver functions for plotting only. % TODO move this somewhere better. 
% yshift_all = [1:nrf]'.* yshift_const; 
transparency = .3 * 150/nrf; % Transparency of X for 100 receiver functions, scaled down if there are more receiver functions. 
if transparency < .15;
    transparency = .15; 
elseif transparency > .5; 
    transparency = .5; 
end
for irf = 1:nrf % Loop over and plot each receiver funcion. 
    yshift = yshift_all(irf); 
    rfi = rf(:,irf);
    hnd_rf = plot(tt, yshift+rfi, 'linewidth', 0.75,...
        'Color', [0,0,0,transparency]); 
end

% ylim([-10*yshift_const, yshift_const * (nrf+10)]); % Find a nice looking y limit. Can play with these values. 
% ylim([-0.25, 0.35]); % Ylimit in amplitude for receiver functions plot. 
lgd = legend([hnd_t_xi(1), hnd_t_0b(1), hnd_t_00(1)], 'Location', 'northeast',...
    'Interpreter','latex', 'fontsize', 12); 

exportgraphics(gcf, sprintf('results/%s/rfs_paper_merged_%s.jpeg',sta_name, sta_name), ...
    'Resolution', res_export); % Save figure. 


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

exportgraphics(gcf, sprintf('results/%s/hk_nice_%s.jpeg',sta_name,sta_name),...
    'Resolution', res_export); 

%% Receiver function plot, with y shift added based on index. 
% Similar to above receiver function plot, but with vertical shift. 

if show_all_plots; % Needs to be updated using the faster plotting code above. 
    figure(202); clf; hold on; 
    set(gcf, 'pos', [2032 336 578 1029]); 
    set(gca, 'LineWidth', 1.5, 'XGrid', 'on', 'XMinorTick', 'on'); box on; 
    xlabel('Time (s)'); 
    title('Receiver functions', 'FontWeight','normal'); 
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
    
    exportgraphics(gcf, sprintf('results/%s/rfs_%s.jpeg',sta_name, sta_name), ...
        'resolution', 250); % save figure. 
end



end 

% save('results/stats/nrf_used.mat', 'nrf_used')