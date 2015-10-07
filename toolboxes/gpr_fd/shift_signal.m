
function odata=shift_signal(idata,dt,nt,omega,eps_mean)
% FFT 2D -> 3D:

T=dt*nt;
Mu0=1.25663706143591730*1e-6; %Magnetic permeability of free space, Vs/Am (4*pi*10^-7 N/A^2)

fftdata=fft2(idata);
fftdata=fftshift(fftdata);

data2=fftdata./sqrt(2*pi*T/(-i*omega*eps_mean*Mu0));

data2=fftshift(data2);
odata=real(ifft2(data2));