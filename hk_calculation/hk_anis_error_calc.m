function [Amap, hVec, kVec, t_pred] = hk_anis_error(...
    rayp, vs, rho, xi, phi, eta, options)  
    arguments
        rayp
        vs
        rho
        xi = 1
        phi = 1 
        eta = 1
        options.phase_wts = [.5, .3, .2]
        options.hBounds = [10, 70]; 
        options.kBounds = [1.5, 2.1]; 
        options.kNum = 201; 
        options.hNum = 200; 
        options.ifplot = false; 
        options.rfinterp = 'cubic'; 
    end

% Handle input values. 
phase_wts = options.phase_wts / sum(options.phase_wts); % Normalize in case these didn't sum to 1 
rayp = rayp/111.1949; 

%% Make grids of H and K
hVec = linspace(options.hBounds(1), options.hBounds(2), options.hNum)';
kVec = linspace(options.kBounds(1), options.kBounds(2), options.kNum) ; 

% Assume ray path change is insignificant when adding anisotropy
incAngVs = hk_rayp2inc(rayp, vs        ) * ones(size(kVec)); 
incAngVp = hk_rayp2inc(rayp, vs * kVec ); %brb timit 8.1e-6

%% Effective vs and vp in anisotropic medium. 
% Solve the Christoffel equation for each possible vp, each k. 
[vs_out] = hk_christof_radial_anis_noang(vs, vs*kVec,... 
            xi, phi, eta, rho, incAngVs); 
vsvAn = vs_out(2,:); 

[vp_out] = hk_christof_radial_anis_noang(vs, vs*kVec,... 
            xi, phi, eta, rho, incAngVp); 
vpAn  = vp_out(1,:); 

%% Calculate phase arrival times. 
% using equations 2-4 from Zhu and Kanamori 2000
kks = sqrt(vsvAn.^-2 - rayp.^2);
kkp = sqrt(vpAn.^-2 - rayp.^2);
t_ps = hVec*(kks - kkp); % using outer product
t_ppps = hVec*(kks + kkp); % using outer product
t_ppss = hVec*(2*kks); % using outer product

% Now that we have timing, go backwards from timing to H using isotropic vs
% and vp. This is how much error there would be if kappa was held constant,
% for various true values of kappa and h. 
vp      = vs*kVec; % Using the isotropic TRUE vp/vs ratio to get vp. 
kks_iso = sqrt(vs.^-2 - rayp.^2);
kkp_iso = sqrt(vp.^-2 - rayp.^2);

hVec_ps = t_ps./(kks_iso - kkp_iso); 
hVec_ppps = t_ppps./(kks_iso + kkp_iso); 
hVec_ppss = t_ppss./(2*kks_iso); 

hVec_iso = + phase_wts(3) .* hVec_ppss ... % use plus for the third phase weight, not minus. Because the timings are negative, it's only the pulse that is negative. 
           + phase_wts(2) .* hVec_ppps ...
           + phase_wts(1) .* hVec_ps; 

fpchange = @(f0, f1)100*(f1-f0)./f0; 

p_change = fpchange(hVec_iso, hVec); 
p_change_t_ps = fpchange(hVec_ps, hVec); 
p_change_t_ppps = fpchange(hVec_ppps, hVec); 
p_change_t_ppss = fpchange(hVec_ppss, hVec); 

figure(1); clf; hold on; 
plot(hVec, p_change); 
xlabel('H true (km)'); 
ylabel('Percent change'); 
title('Percent change H, for different K'); 

figure(2); clf; hold on; 
xlabel('$\kappa$', 'Interpreter', 'latex'); 
ylabel('Percent error H'); 
title('H error analytical', 'FontWeight','normal'); 
LW = 1.5; 
box on; grid on; 
set(gca, 'LineWidth', 1.5); 
plot(kVec, mean(p_change, 1)       , 'DisplayName', 'Average', ...
    'LineWidth', LW); 
plot(kVec, mean(p_change_t_ps,1 )  , 'DisplayName', 'Ps',      ...
    'LineWidth', LW);   
plot(kVec, mean(p_change_t_ppss,1 ), 'DisplayName', 'PpSs' ,   ...
    'LineWidth', LW);   
plot(kVec, mean(p_change_t_ppps,1 ), 'DisplayName', 'PpPs' ,   ...
    'LineWidth', LW);  
legend('location', 'best'); 

end