% sippi_forward_reflection_nmo : AVO reflection seismic forward
%
% Call :
%  [d,forward,prior,data]=sippi_forward_reflection_nmo(m,forward,prior,data)
%
%  m{1}: Vp         % prior{1}.name='vp';
%  m{2}: Vs%        % prior{2}.name='vs';
%  m{3}: Density    % prior{3}.name='rho';
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
%     forward.type='akirichards'; Aki and Richard linear weak contract approximation
%     forward.type='buland_omre'; Aki and Richard linear weak contract
%                                 approximation using Buland and Omre
%                                 formulation
%
%
%
%  % forward.i_trace : number of traces to use. default is to use all traces:
%                    forward.i_trace=1:size(m{1}),2)
%
%
%  See also reflection_convolution_angle
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
    forward.type='akirichard';
    forward.type='zoeppritz';
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

% How many traces should should data be compute for.
if ~isfield(forward,'i_traces')
  forward.i_traces=1:size(m{1},2);
end

%% OBTAIN THE ELASTIC PARAMETERS
% convert m{..} into 'vp', 'vs', 'rho'
if strcmp(lower(prior{1}.name),'vp');
  
  % Vp,Vs,Rho independent prior
    for ip=1:length(prior)
        if strcmp(lower(prior{ip}.name),'vp');
            i_vp=ip;
            vp=m{i_vp}(:,forward.i_traces);
        elseif strcmp(lower(prior{ip}.name),'rho');
            i_rho=ip;
            rho=m{i_rho}(:,forward.i_traces);
        elseif strcmp(lower(prior{ip}.name),'vs');
            i_vs=ip;
            vs=m{i_vs}(:,forward.i_traces);
        end
    end
    
elseif strcmp(lower(prior{ip}.name),'vpvsrho');
    % Cholesky type prior / split vp,vs,rho from 1x[3*nm] vector
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

    
if (strcmp(lower(forward.type),'buland_omre'));
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
    
    if ~isfield(forward,'omre');
      [A,D,W]=buland_omre_setup(forward.vp0,forward.vs0,forward.rho0,ns,forward.angle,forward.wl);
      forward.omre.A=A;
      forward.omre.D=D;
      forward.omre.W=W;
      forward.omre.G=W*A*D;
    end
    
    
    for it=1:ntraces;
      m_all=[vp(:,it);vs(:,it);rho(:,it)];
      d{it}=forward.omre.G*m_all;
    end
    
elseif strcmp(lower(forward.type),'elastic_waveform_2d')
    %% MPM: 
    sippi_verbose(sprintf('%s: type ''%s'' not yet implemented',mfilename,forward.type));
else 
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
    
%else
%  sippi_verbose(sprintf('%s: type ''%s'' not yet implemented',mfilename,forward.type));
end


