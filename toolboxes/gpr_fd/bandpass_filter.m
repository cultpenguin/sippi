

function data_out=bandpass_filter(data,F_stop1,F_pass1,F_pass2,F_stop2,Fs,mode)

% Call: data_out=bandpass_filter(data,F_stop1,F_pass1,F_pass2,F_stop2,Fs,mode);
% * mode: (1): Bandpass, (2): Highpass, (3): Lowpass
% * Fs = 1 / delta_time 

[a b]=size(data);
if a>b
    data=data';
end

A_stop1 = 100;		% Attenuation in the first stopband = 60 dB
%F_stop1 = 1*10^6;%5*10^6;  % Edge of the stopband = 8400 Hz
%F_pass1 = 2*10^6; %10*10^6;	% Edge of the passband = 10800 Hz
%F_pass2 = 150*10^6;	% Closing edge of the passband = 15600 Hz
%F_stop2 = 200*10^6;	% Edge of the second stopband = 18000 Hz
A_stop2 = 100;		% Attenuation in the second stopband = 60 dB
A_pass = 1;		% Amount of ripple allowed in the passband = 1 dB
%Fs = 1/(dt*10^-9);
if mode==1
    BandPassSpecObj=fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2',F_stop1,F_pass1,F_pass2,F_stop2,A_stop1,A_pass,A_stop2,Fs);
    filt = design(BandPassSpecObj,'cheby2');
elseif mode==2
    HighPassSpecObj=fdesign.highpass(F_stop1,F_pass1,A_stop1,A_pass,Fs);
    filt = design(HighPassSpecObj,'cheby2');
elseif mode==3
    LowPassSpecObj=fdesign.lowpass('Fp,Fst,Ap,Ast',F_pass2,F_stop2,A_pass,A_stop2,Fs);
    filt = design(LowPassSpecObj,'cheby2');
else
    disp('Choose another mode')
end

x2=data;
y=filter(filt,x2);
y=filter(filt,fliplr(y));
data_out=fliplr(y);