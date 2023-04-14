% brb2023.04.04 Some examples of the relationship of xi to percent
% anisotropy. Different researchers use different definitions of percent
% anisotropy, so it can be difficult to compare their results. 

%% Lynner et al 2018 used RA = ((Vsh/Vsv)-1)×100).
RA = linspace(-10, 10); 
xi = (RA/100 + 1).^2; 
figure(1); clf; hold on; 
plot(RA, xi); 
xlabel('RA = ((Vsh/Vsv)-1)×100)');
ylabel('\xi'); 


%% Dreiling et al. (2017) used: RA =(𝜉 − 1)×100%
RA = linspace(-15, 15); 
xi = (RA/100 + 1); 
figure(2); clf; hold on; 
plot(RA, xi); 
xlabel('RA = ((Vsh/Vsv)-1)×100)');
ylabel('\xi'); 

%% Dalton et al., 2013: radial anisotropy [2(VSH − VSV)/(VSH + VSV)]

RA = linspace(-.2, .2); 
xi = ( -(RA+2)./(RA-2) ).^(2); 
figure(2); clf; hold on; 
plot(RA*100, xi); 
xlabel('RA = 2 (vsh-vsv)/(vsh+vsv)');
ylabel('\xi'); 
grid on; 