% hk_anis()
nrf = 600; 
rf = randn(nrf,1); 
tt = linspace(0, 40, nrf)'; 
rayp = 6.14; 
vs = 3.3; 
rho = 3.2; 
xi = 1.1; 
phi = 0.9; 
eta = 1; 
hk_anis(rf, tt, rayp, vs, rho, xi, phi, eta); 