% reflectivity_convolution: convolve a reflectivity series with a wavelet
%
% Call: 
%    [seis_data,bigG,G]=reflectivity_convolution(rpp,wavelet,use_method,bigG,G)
% 
%    % INPUT
%    rpp [nt,1]: reflecitivty series
%    wavelet: 
%      wavelet(1).t: time of first wavelet
%      wavelet(1).data: first wavelet
%      wavelet(2).t: time of second wavelet
%      ....
%
%    use_method=1; % convolution 
%    use_method=2; % FFT based convolution
%    use_method=3; % matrix convolution (the convolution matrix is output
%                  % as bigG
%
%    % OUTPUT
%    seis_data: reflection seimic data
%
% See also: reflection_convolution_angle 
%
function [seis_data,bigG,G]=reflectivity_convolution(rpp,wavelet,use_method,bigG,G)

if nargin<3
    use_method=1; % convolution, time domain
    %use_method=2; % convolution, fft domain
    %use_method=3; % matrix convolution
end

if nargin<4,bigG=[];end
if nargin<5,G{1}=[];end

[ns]=size(rpp,1);
nangle=length(wavelet);
seis=zeros(ns,nangle);

sippi_verbose(sprintf('%s: use_method=%d',mfilename,use_method),10);
% loop over angles
for i=1:nangle
    % decide which method to use
    
    if use_method==1,
        % convolution time domains
        for ia=1:nangle
            if length(wavelet)==1
                wl=wavelet(1).data;
                t=wavelet(1).t;
            else
                wl=wavelet(ia).data;
                t=wavelet(ia).t;
            end
            seis_data(:,ia)=real(conv_wl(rpp(:,ia)  ,wl,t));
        end
    elseif use_method==2,
        % fft base convolution
    
    elseif use_method==3,
        %% matrix convolution using fft
        
        %% 1. setup G with wavelet
        if isempty(bigG)
            % compute bigG
            for ia=1:nangle;
                G{ia}=zeros(ns,ns);
                if length(wavelet)==1
                    G{ia}=setup_wavelet_matrix(wavelet(1).data,wavelet(1).t,size(rpp,1));
                else
                    G{ia}=setup_wavelet_matrix(wavelet(ia).data,wavelet(ia).t,size(rpp,1));
                end
                % combine wavelet matrix for each angle into one big matrix
                if ia==1;
                    bigG=G{ia};
                else
                    bigG=blkdiag(bigG,G{ia});
                end
            end
        end
        
        % CONVOLVE
        seis_trace=bigG*rpp(:);
        seis_data=reshape(seis_trace,size(rpp,1),size(rpp,2));
    end
end



