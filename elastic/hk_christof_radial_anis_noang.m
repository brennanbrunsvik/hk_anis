function [velocity] = hk_christof_radial_anis(...
    vs, vp, xi, phi, eta, rho, angAll);
% brb2022.02.25
% This version of code uses explicit expressions for vsv, vsh, and vp from
% the Mavko Rock Physics Handbook. It does not solve for ray orientations,
% which would take much more time. 
% vs, vp: Isotropic portion of velocity
% xi    : radial s anisotropy
% phi   : radial p anisotropy
% eta   : the other anisotropic parameter
% rho   : density
% ang   : ray angle from symmetry axis
% dip_ax: TODO include the dip angle of the symmetry axis, so we can
% account for plunging anisotropy. Should be easy. 
% -------------------------------------------------------------------------
% velocity(i,j): velocity of each quasi wave. 
% i is wave type and j is the ray. Sorted as [vp, vsv, vsh]. 

% get vsv, vsh, vpv, vph, for rays travlling parallel/perpendicular to symmetry direction
[ ys,xs ] = hk_VsvVsh_from_VsXi( vs, xi );
[ yp,xp ] = hk_VpvVph_from_VpPhi( vp, phi );

% Maupin and Park, 2007
A = rho .* xp.^2; % if no anis, xp just P-wave modulus
C = rho .* yp.^2; % if no anis, yp just P-wave modulus
L = rho .* ys.^2; % if no anis, ys just shear modulus
N = rho .* xs.^2; % if no anis, xs just shear modulus
F = eta.*(A - 2*L);% if no anis, is just lambda

% Prep for Voigt notation
c11=A; c12 = A-2*N; c33=C; c44=L; c66=N; c13=F; % Entries to radial anisotropy tensor
c22=c11; c55=c44; c23=c13; % Symmetric entries. 

% % % % Symmetry is in third direction. Would fill matrix as follows. 
% % % CCor = [c11 c12 c13  0  0  0  ; 
% % %         c12 c22 c23  0  0  0  ; 
% % %         c13 c23 c33  0  0  0  ;
% % %          0   0   0  c44 0  0  ; 
% % %          0   0   0   0 c55 0  ; 
% % %          0   0   0   0  0 c66];   

%%% Vp, Vsv, Vsh. for transverse isotropy 
% Eq 2.24-2.2.27 from Mavko Rock Physics Handbook, 3rd edition
theta = angAll; % If symmetry axis is vertical, angAll is the angle between ray propogation and symmetry axis. This would be very easy to change to introduce dipping symmetry axis! 

M   = (   (c11-c44).*sind(theta).^2                                ...
        - (c33-c44).*cosd(theta).^2    ).^2                        ...
        + (c13+c44).^2 .* sind(2*theta).^2                         ; 

Vp  = (   c11.*sind(theta).^2 + c33.*cosd(theta).^2                ...
        + c44 + sqrt(M)                             ).^(1/2)       ...
       .* (2*rho).^(-0.5)                                          ; 

Vsv = (   c11.*sind(theta).^2 + c33.*cosd(theta).^2                ...
        + c44 - sqrt(M)                             ).^(1/2)       ...
       .* (2*rho).^(-0.5)                                          ; 

Vsh = (  (c66.*sind(theta).^2 + c44.*cosd(theta).^2)./rho  ).^(.5) ; 

velocity = [Vp; Vsv; Vsh]; % Put into vector for easy output. 
 
end






















