
function f=disc_freq(N,dt)

% Discrete frequencies:
% Call: f=disc_freq(N,dt);
% N is the number of elements contained in the signal
% dt is the temporal or spatial sampling interval
%
% Knud S. Cordua, 2010
%
% Verification by symetri:
% 
% dt=0.01;
% N=2^nextpow2(1000)+1;
% Tp=1;
% w=gausswavelet(Tp,N,dt);
% wf=fft(w);
% wf=fftshift(wf);
% pwr=real(wf).^2+imag(wf).^2;
% f=disc_freq(N,dt);
% figure, plot(f,pwr)

if ~iseven(1:N)
    % Sampling frequency;
    Fs = 1/dt; 
    f_pos = Fs/2*linspace(0,1,(N-1)/2+1);
    f_neg = Fs/2*linspace(0,1,(N-1)/2+1);
    f=[-fliplr(f_neg(2:end)),f_pos]; % Frequencies
else
    % Sampling frequency;
    Fs = 1/dt; 
    f_pos = Fs/2*linspace(0,1,N/2);
    f_neg = Fs/2*linspace(0,1,N/2+1);
    f=[-fliplr(f_neg(2:end)),f_pos]; % Frequencies
end