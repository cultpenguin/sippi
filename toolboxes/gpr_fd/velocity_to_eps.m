% velocity_to_eps
%
% eps_r is the relative dieletric permittivity
% v is the velocity of the phase (m/ns)
%
% sig is the eletrical conductivity measured in mS/m
% f is the frequency. If f is set to 0 a high frequency approxiamtion is
% applied
%
% (C) Knud Cordua, 2016, Thomas Mejer Hansen, 2016
%

function eps_r=velocity_to_eps(v,sig,f)

if nargin<3
    f=0;
end

% Set the physical constants:
MU0=1.25663706143591730*1e-6; %Magnetic permeability of free space, Vs/Am (4*pi*10^-7 N/A^2)
EPS0=8.85418781762039080*1e-12; %in As/Vm
c0=2.99792458*10^8;
%*-*********-**************************++
velo=v*10^9;
omega=2*pi*f;
mu=MU0;

if f==0
    eps=1./(mu.*velo.^2);
else
    %velo=c0/sqrt(0.5*eps_r*mu_r*(sqrt(1+(sig/(omega*eps))^2)+1));
    sig=sig*10^-3;
    alfa=(sig.*mu.*velo)./2;
    eps=(1./mu).*(1./velo.^2-(alfa./omega).^2);
end

eps_r=eps./EPS0;