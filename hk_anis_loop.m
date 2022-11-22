function [Amap, hVec, kVec, t_pred] = hk_anis_loop(...
    rf, tt, rayp, vs, rho, xi, phi, eta, options)  
    arguments
        rf 
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
    end


nrf = size(rf, 2); 
AmapAll = zeros(nrf,    options.hNum, options.kNum); 
t_pred  = zeros(nrf, 3, options.hNum, options.kNum); 

nrf = length(rayp); 
for irf = 1:nrf; 
    [AmapAll(irf,:,:), hVec, kVec, t_pred(irf,:, :,:)] = hk_anis(...
        rf(:,irf), tt, rayp(1,irf), vs, rho, xi, phi, eta, ...
        'phase_wts', options.phase_wts, ...
        'hBounds', options.hBounds,...
        'kBounds', options.kBounds,...
        'kNum',options.kNum,...
        'hNum',options.hNum, ...
        'ifplot',options.ifplot);  
end

Amap = zeros(size(AmapAll,2), size(AmapAll,3) ) ; 
Amap(:,:) = sum(AmapAll, 1); 
Amap = Amap / max(max(Amap)); 

% interpn(H, K, t_pred00() )

end