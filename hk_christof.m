function [vel, pvecs]=hk_christof(rvec, cc, den)
% This function moved from /Users/brennanbrunsvik/Documents/UCSB/ENAM/VsAnisTomog/synth_split/christof.m

% function [vel, pvecs]=christof(rvec, cc, den);
%
%  Determines velocities and polarizations for wave traveling
%  along ray defined by "rvec", for anisotropic medium with elasticity
%  tensor "cc" and density "den", via Christoffel equation solutions.
%  All coordinates are "local", relative to cc (e.g., (1) = fast axis)
%
% INPUT
%   rvec(): if length(rvec)=3, rvec = cartesian unit vector ray direction
%           if length(rvec)=2, rvec(1) = angle from (1) axis in degrees;
%                             rvec(2) = azimuth from (2)-axis CW in (2,3) plane
%   cc(3,3,3,3) elasticity tensor; same coordinates as rvec  (GPa)
%   den:  density (g/cc)
%
% OUTPUT
%   vel(3): 3 wave velocities; vel(1)=slow-S, vel(2)=fast-S; vel(3)="P" (km/s) 
%  pvecs(3,3):  polarization vectors, each column for each corresponding vel 
%                e.g., pvecs(:,2) is for fast S direction  traveling at vel(2)
%  GAA 7/08

rpd=0.017453292519943;

if (length(rvec)==2)
    th=rvec(1);
    az=rvec(2);
    n=zeros(3,1);
    n(1)=cos(th.*rpd);
    nx=sin(th.*rpd);
    n(3)=nx.*cos(az.*rpd);
    n(2)=nx.*sin(az.*rpd);
elseif (length(rvec)==3)
    n=reshape(rvec,3,1);
    nlen=n'*n;
    if (abs(nlen-1)>.01)
        disp(sprintf('CHRISTOF WARNING: rvec not unit length but %.3f, resizing', nlen));
        n=n./sqrt(nlen);
    end
end
    
% build Christoffel matrix from cc and n
chris=zeros(3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                chris(i,j)=chris(i,j)+cc(k,i,j,l).*n(k).*n(l);
            end
        end
    end
end
chris=chris./den;
[v,d] = eig(chris, 'vector');
[vel,idx]=sort(sqrt(d));
pvecs=v(:,idx);

return