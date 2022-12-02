function [Amap, hVec, kVec, t_pred] = hk_anis(...
    RF, tt, rayp, vs, rho, xi, phi, eta, options)  
    arguments
        RF 
        tt
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
vsvAn = zeros(1, options.kNum); % Row vector. h is collumn vector. 
vpAn  = zeros(1, options.kNum); 

if ~ ((xi==1) && (phi==1) && (eta==1)) ; % Only use time consuming Christoffel equations if there is anisotropy. 
    for ik = 1:options.kNum;
        % Velocities for ray with vs incidence angle
        [velVs,~] = hk_christof_radial_anis(vs, vs*kVec(ik),... 
            xi, phi, eta, rho, incAngVs(ik)); 
        vsvAn(ik) = velVs(2); % Vsv for s ray 
    
        % Velocities for ray with vp incidence angle
        [velVp,~] = hk_christof_radial_anis(vs, vs*kVec(ik),...
            xi, phi, eta, rho, incAngVp(ik)); 
        vpAn (ik) = velVp(1); % Vp for p ray
    end
else % If there's no anisotropy, we already know vp and vsv. 
    vsvAn(:) = vs;
    vpAn (:) = vs * kVec; 
end

%% Calculate phase arrival times. 
% using equations 2-4 from Zhu and Kanamori 2000
kks = sqrt(vsvAn.^-2 - rayp.^2);
kkp = sqrt(vpAn.^-2 - rayp.^2);
t_ps = hVec*(kks - kkp); % using outer product
t_ppps = hVec*(kks + kkp); % using outer product
t_ppss = hVec*(2*kks); % using outer product

% Stored predicted times corresponding to any h, k combo. 
t_pred = zeros(3, size(t_ps,1), size(t_ps,2) ); 
t_pred(1,:,:) = t_ps; t_pred(2,:,:) = t_ppps; t_pred(3,:,:) = t_ppss; 

%% sum weighted contributions from each phase type
% Amap =  phase_wts(1).*interp1(tt,RF,t_ps) ...
%       + phase_wts(2).*interp1(tt,RF,t_ppps) ...
%       - phase_wts(3).*interp1(tt,RF,t_ppss,[],0); % negative phase!

Amap =  phase_wts(1).*interp1(tt,RF,t_ps, options.rfinterp) ...
      + phase_wts(2).*interp1(tt,RF,t_ppps, options.rfinterp) ...
      - phase_wts(3).*interp1(tt,RF,t_ppss,options.rfinterp,0); % negative phase!

end