function [velocity, polarization] = hk_christof_radial_anis(...
    vs, vp, xi, phi, eta, rho, angAll);
% brb2022.02.25
% vs, vp: Isotropic portion of velocity
% xi    : radial s anisotropy
% phi   : radial p anisotropy
% eta   : the other anisotropic parameter
% rho   : density
% ang   : ray angle from symmetry axis
% -------------------------------------------------------------------------
% velocity    : velocity of each quasi wave. Sorted as [vp, vsv, vsh]. 
% polarization: vectors of particle motion polarization. 
%    polarization(:,1) = vector for vp, and so on. 
%    polarization(x, y, z : vp, vsv, vsh)
%    Assume propogation in y and z dimensions only. 
%    Symmetry axis is in positive z direction. 
    
% get vsv, vsh, vpv, vph, for rays travlling parallel/perpendicular to symmetry direction
[ ys,xs ] = VsvVsh_from_VsXi( vs, xi );
% fprintf('\nyx = %1.8f: xs = %1.8f\n', ys, xs)
[ yp,xp ] = VpvVph_from_VpPhi( vp, phi );

% Maupin and Park, 2007
A = rho .* xp.^2; % if no anis, xp just P-wave modulus
C = rho .* yp.^2; % if no anis, yp just P-wave modulus
L = rho .* ys.^2; % if no anis, ys just shear modulus
N = rho .* xs.^2; % if no anis, xs just shear modulus
F = eta.*(A - 2*L);% if no anis, is just lambda

% Insert to Voigt notation matrix
c11=A; c12 = A-2*N; c33=C; c44=L; c66=N; c13=F; % Entries to radial anisotropy tensor
c22=c11; c55=c44; c23=c13; % Symmetric entries. 

% Symmetry is in third direction. 
CCor = [c11 c12 c13  0  0  0  ; 
        c12 c22 c23  0  0  0  ; 
        c13 c23 c33  0  0  0  ;
         0   0   0  c44 0  0  ; 
         0   0   0   0 c55 0  ; 
         0   0   0   0  0 c66];       

% Build full 3x3x3x3 tensor
cc=zeros(3,3,3,3);  
for i=1:3
    for j=1:3
        cc(i,i,j,j)=CCor(i,j);
    end
end
for i=1:2
    for j=(i+1):3
        k=9-j-i;
        c=CCor(k,k);
        cc(i,j,i,j)=c;
        cc(i,j,j,i)=c;
        cc(j,i,j,i)=c;
        cc(j,i,i,j)=c;     
    end
end

% Temp test
% angAll = [1:90]; % Try all incidence angles. 

velocity = zeros(3,length(angAll));
polarization = zeros(3,3,length(angAll)); 
for iRay=[1:length(angAll)];
    ang = angAll(iRay); 

    zerVec = zeros(length(ang),1); 
    propDir = [ zerVec, sind(ang), -cosd(ang)]'; % ray propogation direct (sign doesn't matter). For some reason I drew propogation vector going down, positive is up. 
    [vel,pvecs]=christof(propDir, cc, rho); % Phase velocity and particle motion vectors for three quasi waves. 
    
    % Sort out which is p, sv, sh. 
    pInd = (abs(propDir' * pvecs)); % p is in same direction as propogation
    pInd = find(pInd == max(pInd)); 

    vsvDir  = [zerVec, cosd(ang),  sind(ang)]'; % For figuring out which index is sv
    svInd = abs(vsvDir' * pvecs); 
    svInd = find(svInd == max(svInd)); 

    vshDir = [ zerVec+1, zerVec, zerVec]'; 
    shInd = abs(vshDir' * pvecs); 
    shInd = find(shInd == max(shInd)); 

    % Make sorted vectors of velocity, polarization
    velocityTemp = [vel(pInd), vel(svInd), vel(shInd)]'; 
    polarizationTemp = [pvecs(:,pInd)'; pvecs(:,svInd)'; pvecs(:,shInd)' ]; 
    
    velocity(:,iRay) = velocityTemp; 
    polarization(:,:,iRay) = polarizationTemp; 

    % TODO should also get group velocity (and directions?)
    % See shearer anisotropy chapter for simple formula to go from phase to group.
    % Group is actual ray speed, which is what is measured in many cases. 
end

% figure(1); clf; hold on; 
% lWidth = 2; 
% set(gcf, 'color','white'); 
% box on; 
% grid on; 
% xlabel('Incidence angle'); 
% ylabel('Velocity'); 
% title(sprintf('Velocity(angle). Xi = %1.2f. VSV=%1.8f. VSH=%1.8f',xi,ys,xs)); 
% plot(angAll, velocity(2,:), 'DisplayName', 'Vsv','LineWidth',lWidth); 
% plot(angAll, velocity(3,:), 'DisplayName', 'Vsh','LineWidth',lWidth); 
% plot(angAll, ones(length(velocity),1).* vs, 'DisplayName', 'V','LineWidth',lWidth); 
% leg=legend(); 
% leg.Location = 'best'; 
 
end






















