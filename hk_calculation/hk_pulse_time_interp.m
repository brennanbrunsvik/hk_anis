function [t_terp] = hk_pulse_timing(H, K, t_pred, hterp, kterp)

nrf = size(t_pred,1); 
nh  = size(t_pred, 3); 
nk  = size(t_pred, 4); 

t_terp = zeros(nrf, 3); 
for irf = 1:nrf
    for ipls = 1:3
        t_terp(irf, ipls) = interpn(...
            H, K, reshape(t_pred(irf,ipls,:,:), nh, nk), hterp, kterp); 
    end
end

end