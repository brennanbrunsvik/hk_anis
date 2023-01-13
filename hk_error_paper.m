clear; 
clc; 
restoredefaultpath; 

addpath("elastic/"); 
addpath("hk_calculation/"); 
addpath('~/Documents/repositories/Base_code/colormaps/colormaps_ander_biguri'); 

%%rayp
% rayp = linspace(4, 9, 100); 
rayp = 5; 
vs = 3.8 .* ones(size(rayp));
% vp = 1.75 .* ones(size(rayp));
rho = 2.6 .* ones(size(rayp));
xi = 1.15; 

hk_anis_error_calc(rayp, vs, rho, xi, 1./xi, 1, 'kNum', 5, 'hNum', 6); 