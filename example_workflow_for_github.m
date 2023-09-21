% Do anisotropic HK stacks for a particular station. 
% Make figures as seen in the main paper.  
% This is a simple example script which could be modified to make your own
% anisotropic HK example scripts. 

clear; 
clc; 
restoredefaultpath; 

addpath("./elastic/"); 
addpath("./hk_calculation/"); 

run_version = 'version_1'; % If you want to keep track of a "version" of the code. 
res_export = 200; % Choose base figure resolution DPI

kNum = 150; % How many K to have in grid search. 
hNum = 151;  % How many h to have in grid search. 
hBounds = [10, 55]; % Min and max values of H to include in grid search. 

acparms = [0.85, 1; 0.45, 40]; % Receiver function correlation parameters. [0.85, 1; 0.45, 40] means receiver function should have correlation of 0.85 in first 1 second, 0.45 in first 40 seconds. 

rm_old_fold = true; if rm_old_fold; disp('Removing previous station results. '); end % By default, you should remove the old figures/results before saving the new ones. 

sta_name_all = ["TA.KMSC"]'; % List of stations to run. For the github example, only data for TA KMSC is included. 

% Loop over stations and apply anisotropic HK stack. 
for ista = 1:length(sta_name_all); 
    sta_name = sta_name_all(ista); 
    fprintf("\n ----- %s ---- nsta = %5.0f / %5.0f \n", sta_name, ista, length(sta_name_all))
    
    if rm_old_fold & (exist(sprintf('./results/%s', sta_name))~=0); 
        rmdir(sprintf('./results/%s', sta_name), 's'); 
    end
    
    sta_name_split = sta_name.split('.'); 
    sta = sta_name_split(2); 
    net = sta_name_split(1); 
    
    %% Prepare model. 
    % Included is a vectors of model parameters with depth. 
    % If you use such a format, you need to break them down to average
    % crustal values. hk_average_crust() does that. If you already have
    % average velocity and xi, skip this step and provide [rhoav, vsav, vpav, xiav, phiav, etaav]

    % Load an example model. In order to modify and use  
    thing = load(sprintf('models/model_%s_%s.mat', lower(net), lower(sta)) ); % Loads variables rf, tt, eta, phi, rayp, etc. 
    eta = thing.eta; % Anisotropy parameter
    phi = thing.phi; % Anisotropy parameter
    rayp = thing.rayp; % Ray parameter
    rf = thing.rf; % Receiver functions
    rho = thing.rho; % Density
    tt = thing.tt; % Time for receiver functions
    vp = thing.vp; 
    vs = thing.vs; 
    xi = thing.xi; % Radial anisotropy
    z = thing.z; % Depth where xi, vp, etc., is known. 
    zmoh = thing.zmoh; % If you want to apply earth flattening transformation, use a preliminary zmoh to decide the depth to which to do this. 

    % Get average crustal values
    [rhoav, vsav, vpav, xiav, phiav, etaav] =  hk_average_crust(...
        rho, vs, vp, xi, phi, eta, z, zmoh);  
    
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
        ac = nan(nrf, nrf); % correlation matrix. 
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
    
    
    if ~ exist(sprintf('./results/%s', sta_name)); 
        mkdir(sprintf('./results/%s', sta_name)); 
    end
    
    subplot(2,1,2); hold on; 
    title('Kept'); 
    plot(tt, rf); 
    xlim([-5, 40]); 
    
    % Sort by ray parameter. Makes plotting more obvious. 
    [rayp, sortp] = sort(rayp); 
    rf = rf(:,sortp); 
    
    %% Make receiver functions. Do isotropic and anisotropic. 
    % This is the most important part of the script. Use debugging to
    % follow what happens in hk_anis_loop. 
    % Exi: include xi in receiver functions. E00: standard process without xi.  
    [Exi, H, K, t_pred_all_xi] = hk_anis_loop(...
        rf, tt, rayp, vsav, rhoav, xiav, phiav, etaav, ...
                'kNum', kNum, 'hNum', hNum, 'normalize', false, ...
                'hBounds', hBounds); % accounting for xi
    [E00, H, K, t_pred_all_00] = hk_anis_loop(...
        rf, tt, rayp, vsav, rhoav,    1,     1,     1, ...
                'kNum', kNum, 'hNum', hNum, 'normalize', false, ...
                'hBounds', hBounds); ; % not accounting for anisotropy
    
    %% Get the "Best" h and kappa values. Calculate pulse timings. 
    [ihmax, ikmax] = find(Exi == max(Exi ,[], 'all')); % H and K that give highest energy when anisotropy is accounted for. 
    kmax = K(ikmax); 
    hmax = H(ihmax); 
    
    [ihmaxiso, ikmaxiso] = find(E00 == max(E00 ,[], 'all')); % H and K that give highest energy when anisotropy is ignored.
    kmaxiso = K(ikmaxiso); 
    hmaxiso = H(ihmaxiso); 
    
    Exi_max = max(Exi, [], 'all'); % Highest energy in anisotropic and anisotropy stacks. 
    E00_max = max(E00, [], 'all'); 
    
    % Calculate pulse timings. 
    t_pred_xi = hk_pulse_time_interp(H, K, t_pred_all_xi, hmax, kmax); % using xi
    t_pred_00 = hk_pulse_time_interp(H, K, t_pred_all_00, hmax, kmax); % O anisotropy, 
    t_pred_0b = hk_pulse_time_interp(H, K, t_pred_all_00, hmaxiso, kmaxiso); % 0 anisotropy but best model
    
    fprintf('\nHmax_xi = %1.3f, Hmax_iso = %1.3f. dH = %1.3f. percent dH = %1.3f\n', ...
        hmax, hmaxiso, hmax-hmaxiso, (hmax-hmaxiso)/hmaxiso*100)
    fprintf('\nKmax_xi = %1.3f, Kmax_iso = %1.3f. dK = %1.3f. percent dK = %1.3f\n', ...
        kmax, kmaxiso, kmax-kmaxiso, (kmax-kmaxiso)/kmaxiso*100)
    fprintf('\nRatio of E_xi over E_iso = %1.2f\n', Exi_max/E00_max)
    
    save(sprintf('results/%s/hk_summary.mat', sta_name), 'hmax', 'hmaxiso', 'kmax', 'kmaxiso', 'Exi_max', 'E00_max', 'n_rf', 'run_version', 'vsav', 'rhoav', 'xiav', 'phiav', 'etaav'); 
    
    %% All code below is to make figures. 
    rf = rf * 2.7; % Scale receiver functions, for plotting. 
    
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
    nxt = 5; % Number of x tiles for tiledlayout
    nyt = 10; % y tiles for tiledlayout
    xy_to_t = @(x,y)(y-1)*nyt+x; % A hack so we can fill multiple tiles with one figure. 
    tiledlayout(nxt, nyt, 'TileSpacing','compact');
    ax1 = nexttile(xy_to_t(1,1), [5,6]); hold on; 
    ax2 = nexttile(xy_to_t(7,1), [5,4]); hold on; 
    
    each_ax = [ax1, ax2]; % Make each axis look pretty. 
    for ieach_ax = 1:length(each_ax); 
        axes(each_ax(ieach_ax)); 
        set(gca, 'LineWidth', 1.5); box on; 
    end
    
    axes(ax1); 
    xlabel('Time (s)'); 
    ylabel('Ray parameter (km/s)'); 
    title('Receiver functions', 'FontWeight','normal', 'Interpreter', 'latex');  
    text(xlbl, ylbl, '(a)', 'Units','normalized'); 
    
    axes(ax2); 
    xlabel('\kappa'); 
    ylabel('H (km)'); 
    grid on;
    title('$H-\kappa$ stack', 'fontweight', 'normal', 'Interpreter', 'latex'); 
    text(xlbl, ylbl, '(b)', 'Units','normalized'); 
    
    % HK stack 
    axes(ax2); 
    fhand_norm = @(inval)inval ./ max(max(inval)); % Return normalized inval 
    set(gca,'ydir', 'reverse');
    plt_ylim = hmax + 5.*[-1, 1]; % Shift limits of plots to within range of H and K. 
    plt_xlim = kmax + .25.*[-1, 1]; 
    
    % If plot has a fairly high/low k or h max, then shift the center of the plot toward the k or h. 
    if max(plt_xlim) > max(K); 
        plt_xlim = plt_xlim - (max(plt_xlim) - max(K)); % Make sure limits aren't out of bounds
    end
    ylim(plt_ylim); 
    xlim(plt_xlim); 
    
    lvl_cnt = [0.8, 0.95]; % Contour percentages for HK stack plot. Level is percent compared to maximum of either HK stack. 
    LW = 1; % Linewidth. 
    
    LW_scale = 1.75; % Increase linewidth by this much for contour lines
    [~,hnd_xistart] = contour(K, H,  Exi,...
        lvl_cnt * max(Exi, [], 'all'), 'LineWidth', LW*LW_scale, 'color', c_t_with_xi.*cnt_mod, ...
        'DisplayName', sprintf('\\xi = %1.2f', xiav ) ); 
    [~,hnd_xiend  ] = contour(K, H,  E00,...
        lvl_cnt * max(E00, [], 'all'), 'LineWidth', LW*LW_scale, 'color', c_t_iso.*cnt_mod, ...
        'DisplayName', sprintf('\\xi = %1.2f', 1.00 )); 
     
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
    
    t_plot_xi_all = zeros(nrf, 3); 
    t_plot_00_all = zeros(nrf, 3); 
    t_plot_0b_all = zeros(nrf, 3); 
    y_plt_xi_all  = zeros(nrf, 3); 
    y_plt_00_all  = zeros(nrf, 3); 
    y_plt_0b_all  = zeros(nrf, 3); 
    for irf = 1:nrf % Loop over each receiver function and plot the timing of the Ps etc. pulses. 
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
    
    transparency = .3 * 150/nrf; % Transparency of X for 100 receiver functions, scaled down if there are more receiver functions. 
    if transparency < .15;
        transparency = .15; 
    elseif transparency > .5; 
        transparency = .5; 
    end
    
    % Loop over and plot each receiver funcion. 
    for irf = 1:nrf 
        yshift = yshift_all(irf); 
        rfi = rf(:,irf);
        hnd_rf = plot(tt, yshift+rfi, 'linewidth', 0.75,...
            'Color', [0,0,0,transparency]); 
    end
    
    [hLg, icons]=legend([hnd_t_xi(1), hnd_t_0b(1), hnd_t_00(1)], 'Location', 'northeast',...
        'Interpreter','latex', 'fontsize', 12); 
    icons = findobj(icons,'Type','patch'); % Increase scatter sizes in the legend.  
    icons = findobj(icons,'Marker','none','-xor'); % Increase scatter sizes in the legend. 
    set(icons(1:3),'MarkerSize',30); % Increase scatter sizes in the legend. 
    
    set(gcf, 'Renderer', 'painters'); 
    exportgraphics(gcf, sprintf('results/%s/rfs_paper_merged_%s.jpeg',sta_name, sta_name), ...
        'Resolution', res_export); % Save figure. 

end 