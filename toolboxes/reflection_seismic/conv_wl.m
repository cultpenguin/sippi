% conv_wl: Convolute wavelet with time seriesconvolution model
%
%  seis=conv_wl(rpp,wl,wl_time,useFFT);

% perhaps a problem for wavelets that are not symmetroix?
function seis=conv_wl(rpp,wl,wl_time,useFFT);
if nargin==2
    wl_time=[0:1:length(wl)];
end
if nargin<4, useFFT=0;end

if useFFT
    seis=conv_fft(rpp,wl);do_Trim=1;
    seis=conv_fft(rpp,wl,'same');do_Trim=1;
else
    seis=conv(rpp,wl);doTrim=1;
    seis=conv(rpp,wl,'same');doTrim=0;
end

if doTrim==1;
    ns=length(rpp);
    % trim the data
    t0=find(wl_time==0);
    if isempty(t0)
        t0=find(min(abs(wl_time))==abs(wl_time));
    end
    seis=seis(t0:1:length(seis));
    seis=seis(1:ns);
end

doPlot=0;
if doPlot==1
    
    ii=1:1:length(rpp);
    subplot(1,3,1);
    plot(rpp,ii,'k-');
    set(gca,'xlim',[-1.1 1.1])
    grid on
    
    subplot(1,3,2);
    plot(wl,wl_time,'g-')
    grid on
    
    subplot(1,3,3);
    plot(seis,1:1:length(seis),'g-')
    grid on;
end



