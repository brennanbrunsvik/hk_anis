function [rhoav, vsav, vpav, xiav, phiav, etaav] = hk_average_crust(...
    rho, vs, vp, xi, phi, eta, z, zmoh)
% If you have a complex crustal model, this make average crustal values
% that are compatible with the Zhu and Kanamori 2000 formulae. 

% Do Earth flattening transformation to slightly improve accuracy. 
% See Peter Shearer Intro Seismology textbook. 
a = 6371; % Earth radius
zflt = @(z)-a .* log((a-z)./a); 
vflt = @(z,v)a./(a-z) .* v; 

z_orig    = z; 
% zmoh_orig = zmoh; 
z         = zflt(z   ); 
zmoh      = zflt(zmoh); 
vs        = vflt(z_orig, vs); 
vp        = vflt(z_orig, vp); 

% Figure out layer spacing and which layers are in the crust. 
isCrust = z <= zmoh; 
isCrustF = find(isCrust); 
if z(isCrustF(end-1)) == z(isCrustF(end)); % Handle possible duplicated moho depth for crust and mantle. e.g. at 45 km, could have crustal vs in one layer and mantle vs in lower layer. 
    isCrustF(end) = []; % Remove mantle value
    isCrust(:) = false; 
    isCrust(isCrustF) = true; 
end
dz = zeros(size(z)); % For weighted mean, integral style. 
dz(1:end-1,1) = dz(1:end-1,1) + 0.5 .* diff(z);
dz(2:end  ,1) = dz(2:end  ,1) + 0.5 .* diff(z); 
dzC = dz(isCrust); 

% Function to calculate mean accounting for variable dz spacing. Simple integral mean. 
meanInt = @(m,dz)sum(dz.*m)/sum(dz); 
% vrAv = @(m,dz)( (1./meanInt(1./m,dz) + meanInt(m,dz) ) / 2 ); % Voigt Reuss Hill average. 
vrAv = @(m,dz)( (1./meanInt(1./m,dz) ) ); % This averaging approach is slightly more consistent with propmat than the voigt reuss. Makes sense... it's all about travel times. 


rhoav  = meanInt(rho(isCrust), dzC); 
etaav  = meanInt(eta(isCrust), dzC); 
xiav   = meanInt(xi (isCrust), dzC); 
phiav  = meanInt(phi(isCrust), dzC); 
vsav   = vrAv(vs(isCrust), dzC); 
vpav   = vrAv(vp(isCrust), dzC); 

end
