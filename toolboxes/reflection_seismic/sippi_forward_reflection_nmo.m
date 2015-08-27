% sippi_forward_reflection_nmo : AVO reflection seismic forward
%
% Call :
%  [d,forward,prior,data]=sippi_forward_reflection_nmo(m,forward,prior,data,id,im)
%
%  m{1}: Vp
%  m{2}: Density
%  m{3}: Vs
%  d : d{i} AVO ANGLE GATHER, for i=1:size(m{1},2)
%
%  % define angle gathers 
%    forward.angle: array of incidence angles
%    forward.freq: attay of the dominant frequency of Ricker wavlet
%    forward.wl: A wavelet for each angle (overrides forward.freq)
%        forward.wl(1).t: time for wavelet 1
%        forward.wl(1).data: wavelet 1
%        forward.wl(2).t: time for wavelet 2
%        forward.wl(2).data: wavelet 2
%        ...
%
%  % forward types
%     % linear convolution
%     forward.type='shuey';           % 3 term Shuey
%     forward.type='shuey_2_term';    % 2 term Shuey
%     forward.type='full_zoeppritz';    % 2 term Shuey
%     forward.type='weak_contrast'; Aki and Richard liear weak contract approximation
%
%     forward.convolution_type='full_zoeppritz';  % Full zoepprit
%
%
%
%
function [d,forward,prior,data]=sippi_forward_reflection_nmo(m,forward,prior,data,id,im)

if nargin<6, im=1;end
if nargin<5, id=1;end
if nargin<4, data{1}.null='';end

if ~isfield(forward,'parameter')
    forward.parameter='time';
    % forward.parameter='depth';
end

if ~isfield(forward,'type')
    forward.type='shuey';
end

if ~isfield(forward,'conv_method')
    forward.conv_method=1; % conv, time
    forward.conv_method=2; % conf fft
    forward.conv_method=3; % matrix conv
end

% do we make use of log-space?
if ~isfield(forward,'log')
    forward.log=0;
end

%% OBTAIN THE ELASTIC PARAMETERS
% convert m{..} into 'vp', 'vs', 'rho'
if strcmp(lower(prior{1}.name),'vp');
    % Vp,Vs,Rho independent prior
    for ip=1:length(prior)
        if strcmp(lower(prior{ip}.name),'vp');
            i_vp=ip;
            vp=m{i_vp};
        elseif strcmp(lower(prior{ip}.name),'rho');
            i_rho=ip;
            rho=m{i_rho};
        elseif strcmp(lower(prior{ip}.name),'vs');
            i_vs=ip;
            vs=m{i_vs};
        end
    end
    
elseif strcmp(lower(prior{ip}.name),'vpvsrho');
    % Cholesky type prior
end


%% CHECK THAT WAVELET IS DEFINED
if ~isfield(forward,'wl');
    % nowavelet is set
    if ~isfield(forward,'freq');
        forward.freq=ones(1,length(forward.angle)).*[35];
    end
    txt=sprintf(' %f',forward.freq);
    disp(sprintf('%s: No wavelet is set.',mfilename));
    disp(sprintf('%s: Using Ricker with center freq=%s',mfilename,txt));
    
    if ~isfield(forward,'t')
        disp(sprintf('%s: please set the time in "forward.t"',mfilename));
    end
    
    % compute Ricker wavelets for each angle
    for i=1:length(forward.freq);
        dt=forward.t(2)-forward.t(1);
        if (length(forward.freq)==1);
            freq=forward.freq(1);
        else
            freq=forward.freq(i);
        end
        [w,t]=rickerwavelet(freq,dt);
        forward.wl(i).data=w;
        forward.wl(i).t=t;
    end
end
%% COMPUTE DATA
% for each forward.type we 
% 1 check for log- and depth-parameterizations
%

if (strcmp(lower(forward.type),'full_zoeppritz')|strfind(lower(forward.type),'shuey'));
    % CHANGE TO LOOP OVER X AND Y TO USE WITH 3D DATA
    if strcmp(forward.parameter,'depth');
        [vp,vs,rho]=depth_to_time(vp,vs,rho,forward.x,forward.y,forward.t);
    end
    if forward.log==1
        seis=reflection_convolution_angle(exp(vp),exp(vs),exp(rho),forward.angle,forward.wl,forward.type);
    else
        seis=reflection_convolution_angle(vp,vs,rho,forward.angle,forward.wl,forward.type);
    end
    
    % output each NMO gather in a data structure
    for ix=1:size(vp,2);
        d{ix}=seis(:,ix);
    end
    
elseif (strcmp(lower(forward.type),'weak_contrast'));
    % Aki and Richard Weak Contrast approximation
    % following Buland and Omre (2003)
    % the backgrond model vp0, vs0, rho0 are set to constant values
    % vs^2/vp^2 is set to a constant
    %
    % d = WADm
    
    % check for log-parameterization
    if forward.log==0
        vp=log(vp);
        vs=log(vs);
        rho=log(rho);
    end
    if ~isfield(forward,'vp0'); forward.vp0=mean(vp(:));end
    if ~isfield(forward,'rho0'); forward.rho0=mean(rho(:));end
    if ~isfield(forward,'vs0'); forward.vs0=mean(vs(:));end
    [ns,ntraces]=size(vp);
    nangle=length(forward.wl);
    if ~isfield(forward,'omre');
        vs2vp2=exp(forward.vs0).^2./exp(forward.vp0).^2;
        forward.omre.a_vp=0.5*(1+tan(pi*forward.angle/180).^2);
        forward.omre.a_vs=-4*(vs2vp2)*sin(pi*forward.angle/180).^2;
        forward.omre.a_rho=.5*(1+forward.omre.a_vs);
        forward.omre.Asmall=[forward.omre.a_vp(:) forward.omre.a_vs(:) forward.omre.a_rho(:)];
        na=length(forward.wl);
        for ia=1:na;
            A1=eye(ns).*forward.omre.Asmall(ia,1);
            A2=eye(ns).*forward.omre.Asmall(ia,2);
            A3=eye(ns).*forward.omre.Asmall(ia,3);
            
            if ia==1
                A=[A1 A2 A3];
            else
                A=[A;A1 A2 A3];
            end
        end
        forward.omre.A=A;
    end
    
    %keyboard
    %zoeppritz_aki_richards_timeseries(vp,vs,rho,forward.angle);
    
    for it=1:ntraces;
        
      m_all=[vp(:,it);vs(:,it);rho(:,it)];
      
      if ~isfield(forward.omre,'D');
        % setup differential operatroe
        n=size(vp,1);
        %D_small=zeros(n,n);for i=1:(n-1);D_small(i,[1:2]+i-1)=[-1 1];end %%pad top
        D_small=zeros(n,n);for i=2:n;D_small(i,[1:2]+i-2)=[-1 1];end % pad bot
        D=blkdiag(D_small,D_small,D_small);
        forward.omre.D=D;
      end
      m_diff=forward.omre.D*m_all;
      % alternatively without ysing difference operator
      %m_diff = [0;diff(vp(:,it));0;diff(vs(:,it));0;diff(rho(:,it))];
      
        
      
      do_forward=1;
      if isfield(forward.omre,'WA')
        % linear forward matrix allready available
        % seis_data = W*A*D*m_all = [W*A]*[D*m_all];
        seis_data=forward.omre.WA*m_diff;
        do_forward=0;
      end
      if do_forward==1
        rpp = reshape(forward.omre.A*m_diff,ns,nangle);
        [seis_data,forward.omre.W,forward.omre.Wsmall]=reflectivity_convolution(rpp,forward.wl,forward.conv_method);
        % store the product of W and A
        forward.omre.WA=forward.omre.W*forward.omre.A;
      end
      
      d{it}=seis_data(:);
    end
    
    
elseif strcmp(lower(forward.type),'elastic_waveform_2d')
    %% MPM
end


