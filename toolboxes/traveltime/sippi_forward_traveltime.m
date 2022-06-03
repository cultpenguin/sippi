% sippi_forward_traveltime Traveltime computation in SIPPI
%
% Call :
%   [d,forward,prior,data]=sippi_forward_traveltime(m,forward,prior,data)
%
%   forward.type determines the method used to compute travel times
%   forward.type='ray_2d';    % raytracing 2D linear forward
%   forward.type='ray';       % ray (obtained using 'eikonal')
%   forward.type='fat';       % fat ray (obtained using 'eikonal')
%   forward.type='eikonal';   % The Eikonal solution 
%   forward.type='born';      % The Born approximation
%   forward.type='fd';        % Full waveform modeling
%
%   forward.sources [ndata,ndim]: Source locations
%   forward.receivers [ndata,ndim]: Receiver locations
%
%   -->[only related to kernels 'ray' and 'fat']
%   forward.linear : [0] a linear kernel is computed, based on the current velocity model
%                    [1] a linear kernel is computed only once, based on
%                    the velocity field defined in forward.linear_m;
%
%
%   forward.is_slowness: [0] Assumes prior on model parameters are defined
%                            in VELOCITY  domain (default)
%   forward.is_slowness: [1] Assumes prior on model parameters are defined
%                            in slowness domain
%
%   forward.freq : [scalar] Signal frequency,m used to define the width of
%                  the kernels forward.G (only kernels based on 'eikonal'
%                  type kernels)
%
%   forward.linear_m: the reference velocity field, for a linear forward
%                     operator (forward.G) will be computed.
%                     Can be eithe a scalar (constant velocity field) or
%                     the same size as the the velcity model 'm'.
%
%   forward.normalize_vertical [1]: Normalize sensitivitykernel by
%                                   in vertical slices
%                              [0]: No normalization
%
%   forward.alpha [1]: alpha value for munk_fresnel_2d, munk_fresnel_3f
%
%   -->[only related to kernels 'ray_2d']
%     forward.r=1; % Refinement factor
%     forward.norm=0; % normalize sum of each row to 1
%     (see also: ray_kernel_2d.m)
%
%
%   See also: munk_fresnel_2d, munk_fresnel_3d, tomography_kernel, eikonal,
%   eikonal_traveltime, ray_kernel_2d
%
function [d,forward,prior,data]=sippi_forward_traveltime(m,forward,prior,data,id,im)


%if nargin<4;    forward.null='';end
if nargin<4;    data{1}.null='';end
if nargin<5;    id=1;end
if nargin<6;
    if isfield(forward,'im')
        im=forward.im;
    else
        im=1;
    end
end

if isfield(data{id},'sources')
    forward.sources=data{id}.sources;
end
if isfield(data{id},'receivers')
    forward.receivers=data{id}.receivers;
end

S=forward.sources;
R=forward.receivers;


if isfield(forward,'x');
    x=forward.x;
else
    x=prior{im}.x;
end
if isfield(forward,'y');
    y=forward.y;
else
    y=prior{im}.y;
end
if isfield(forward,'z');
    z=forward.z;
else
    z=prior{im}.z;
end

if ~isfield(data{id},'i_use');
    %data{id}.i_use=1:size(data{id}.d_obs,1);
    data{id}.i_use=1:size(forward.sources,1);
end

if ~isfield(forward,'normalize_vertical');forward.normalize_vertical=1;;end
if ~isfield(forward,'alpha');forward.alpha=1;;end
if ~isfield(forward,'type');forward.type='eikonal';;end
if ~isfield(forward,'is_slowness');forward.is_slowness=0;;end

% data response
if strcmp(forward.type,'eikonal')
    % EIKONAL
    if ~isfield(forward,'eikonal_type');
        forward.eikonal_type=1; % FM
        %forward.eikonal_type=2; % NFD
    end
    if forward.is_slowness==1;
        % SLOWNESS
        d{id}=eikonal_traveltime(x,y,z,1./m{im},forward.sources,forward.receivers,data{id}.i_use,forward.eikonal_type);
    else
        % VELOCITY
        d{id}=eikonal_traveltime(x,y,z,m{im},forward.sources,forward.receivers,data{id}.i_use,forward.eikonal_type);
    end
elseif (strcmp(forward.type,'ray')|strcmp(forward.type,'fat'));
    % RAY APPROXIMATION   
    if ~isfield(forward,'linear');forward.linear=0;end
    if (~isfield(forward,'G')|forward.linear==0);
        % ONLY COMPUTE KERNEL IF linear=1 OR IF G DOES NOT EXIST
        if ~isfield(forward,'freq'), forward.freq = 0.1;end
        T=1./ forward.freq;
        
        if ~isfield(forward,'linear_m'),
            if ~isfield(forward,'m0'),
                forward.m0=mean(m{im}(:));
            end
            forward.linear_m = m{im}.*0+forward.m0;
        end
        if forward.linear==1
            m_use=forward.linear_m;
        else
            m_use=m{im};
        end
        
        if forward.is_slowness==1;
            [K,RAY,Gk,Gray]=tomography_kernel(1./m_use,x,y,z,S(data{id}.i_use,:),R(data{id}.i_use,:),T,forward.alpha,forward.normalize_vertical);
        else            
            [K,RAY,Gk,Gray]=tomography_kernel(m_use,x,y,z,S(data{id}.i_use,:),R(data{id}.i_use,:),T,forward.alpha,forward.normalize_vertical);
        end
        
        if strcmp(forward.type,'ray')
            forward.G=Gray;
        else
            forward.G=Gk;
        end
    end
    
    if length(data{id}.i_use)==size(forward.G,1)
        ig=1:size(forward.G,1);
    else
        ig=data{id}.i_use;       
    end
    if forward.is_slowness==1
        d{id}=forward.G(ig,:)*m{im}(:);
    else
        s=1./m{im}(:);
        d{id}=forward.G(ig,:)*s;
    end
elseif strcmp(forward.type,'ray_2d');
    % Straight ray approximation using ray_kernel_2d
    
    if ~isfield(forward,'r');
        forward.r=1; % Refinement factor
    end
    if ~isfield(forward,'norm');
        forward.norm=0; % normalize
    end
    
    if (~isfield(forward,'G'));
    
    if ~isfield(prior{im},'ndim');
        prior=sippi_prior_init(prior);
    end
    if prior{im}.ndim~=2
        sippi_verbose(sprintf('%s: ''%s'' only works in 2D',mfilename,forward.type));
        return
    end    
    dx=prior{im}.x(2)-prior{im}.x(1);
    dy=prior{im}.y(2)-prior{im}.y(1);
    if ((dx-dy)./dx)>1e-15
        sippi_verbose(sprintf('%s: ''%s'' only works in case dx=dy (%g=%g)',mfilename,forward.type,dx,dy));
        return
    end
    ant_pos=[forward.sources forward.receivers];
    ant_pos(:,1)=ant_pos(:,1)-prior{im}.x(1)+dx/2;
    ant_pos(:,3)=ant_pos(:,3)-prior{im}.x(1)+dx/2;
    ant_pos(:,2)=ant_pos(:,2)-prior{im}.y(1)+dx/2;
    ant_pos(:,4)=ant_pos(:,4)-prior{im}.y(1)+dx/2;
    forward.G=ray_kernel_2d(ant_pos,length(prior{im}.y),length(prior{im}.x),dx,forward.r,forward.norm);
    end
    
    if length(data{id}.i_use)==size(forward.G,1)
        ig=1:size(forward.G,1);
    else
        ig=data{id}.i_use;       
    end
    if forward.is_slowness==1
        d{id}=forward.G(ig,:)*m{im}(:);
    else
        s=1./m{im}(:);      
        d{id}=forward.G(ig,:)*s;
    end
    
elseif strcmp(forward.type,'born');
    disp('born')
    % BORN APPROXIMATION
    % LIU ET AL (2009) for seismic waves
    % Buursink et al () for electromagnetic waves
    if ~isfield(forward,'linear_m'),
          if ~isfield(forward,'m0'),
              forward.m0=mean(m{im}(:));
          end
        forward.linear_m = m{im}.*0+forward.m0;
    end
    if ~isfield(forward,'linear'),forward.linear=1;end
    if ~isfield(forward,'freq');
        forward.freq=0.1;
        disp(sprintf('''forward.freq'' field is not set !!',mfilename))
        disp(sprintf('''Using forward.freq''=%g field is not set !!',mfilename,forward.freq))
    end
    if forward.linear==1
        m_use=forward.linear_m;
        use_eikonal=0; % fastest born kernel
    else
        m_use=m{im};
        use_eikonal=1; % slower born kernel, but valid for small vel-contrast
    end
    if ~isfield(forward,'G');
        forward.G=zeros(length(data{id}.i_use),prod(prior{im}.dim));
        j=0;
        for i=data{id}.i_use;
            j=j+1;
            progress_txt(i,size(forward.sources,1),'setting up kernel');
            
            if forward.is_slowness==1;
                [kernel,L,L1_all,L2_all]=kernel_buursink_2d(1./m_use,x,y,forward.sources(i,:),forward.receivers(i,:),forward.freq,[],use_eikonal);
            else
                [kernel,L,L1_all,L2_all]=kernel_buursink_2d(m_use,x,y,forward.sources(i,:),forward.receivers(i,:),forward.freq,[],use_eikonal);
            end
            
            
            forward.G(j,:)=kernel(:);
        end
    end
    if forward.is_slowness==1
        d{id}=forward.G*m{im}(:);
    else
        s=1./m{im}(:);
        d{id}=forward.G*s;
    end

elseif strcmp(forward.type,'fd');
    %% ERNST FW CODE WITH KNUDS PARRALLEL MATLAB INTERFACE
    
    %% full waveform modeling
    if ~isfield(forward,'fd'); forward.fd.null='';end
    
    forward.fd.output_type='trace'; % 'trace' or 'gather'
    forward.fd.ant_pos=[forward.sources forward.receivers];
    %forward.fd.ant_pos=forward.fd.ant_pos(1:1:351,:);
    
    forward.fd.output_it=1; % output every 'output_ti' samples
    if ~isfield(forward.fd,'t');forward.fd.t=1e-7; end;% Frequency
    if ~isfield(forward.fd,'sig');forward.fd.sig=3; end;% constant
    if ~isfield(forward.fd,'debug');forward.fd.addpar.debug=-1;end
    %if ~isfield(forward.fd,'sig');forward.fd.addpar.cores=4;end
    %forward.fd.addpar.Tg=400*10^6
    
    
    
    
    %% convert from velocity to dielectric permittivity
    forward.m_fd{1}=velocity_to_eps(m{1});
    % FD forward response
    [d_fd,forward.fd]=sippi_forward_gpr_fd(forward.m_fd,forward.fd,prior); 
    forward.fd.d=d_fd;
    %% first arrival picking    
    if ~isfield(forward,'fa'); forward.fa.null='';end
    if ~isfield(forward.fa,'ref_to');forward.fa.ref_t0=0;end
    if ~isfield(forward.fa,'doPlot');forward.fa.doPlot=0;end
    if ~isfield(forward.fa,'use_method');forward.fa.use_method=2;end
    forward.fa.wf_time=forward.fd.time;
    if ~isfield(forward.fa,'ref_trace');
        forward.fa.ref_trace=forward.fd.d{1};
    end
    
    forward.fa.doPlot=0;    
    for i=1:size(forward.fd.ant_pos,1);
        %try
        wf_data=forward.fd.d{i};
        [tt_pick(i),t(i)]=pick_first_arrival(wf_data,forward.fa.ref_trace,forward.fa.ref_t0,forward.fa.doPlot,forward.fa.wf_time,forward.fa.use_method);
        %catch
        %    i
        %    keyboard
        %end
          
    end
    t=t*1e+9;
    
    %% optionally compute time shift to apply to first arrival pick 
    if ~isfield(forward.fa,'t_shift');
        sippi_verbose(sprintf('%s: computing reference t0 for traveltime picking',mfilename))
        % compute FD in one data set in a homo
        v0=mean(m{1}(:));
        eps0=velocity_to_eps(v0);
        forward.m_fd_0{1}=forward.m_fd{1}.*0+eps0;
        
        f=forward.fd;
        %f.ant_pos=forward.fd.ant_pos(1,:);
        pause(.5); %sometime necessary
                
        [fd_0]=sippi_forward_gpr_fd(forward.m_fd_0,f,prior);
        forward.fd.d_fd_0=fd_0;
        for i=1:size(forward.fd.ant_pos,1);
            try 
                wf_data=forward.fd.d_fd_0{i};
                [tt_pick_0(i),t_0(i)]=pick_first_arrival(wf_data,forward.fa.ref_trace,forward.fa.ref_t0,forward.fa.doPlot,forward.fa.wf_time,forward.fa.use_method);
            catch
                t_0(i)=0;
                sippi_verbose(sprintf('%s: something went wrong',mfilename))
                %keyboard
            end
        end
        t_0=t_0*1e+9;
    
        % compute traveltimes analytically       
        for i=1:size(f.ant_pos,1);
            t_ex(i) = sqrt((f.ant_pos(i,3)-f.ant_pos(i,1))^2+(f.ant_pos(i,4)-f.ant_pos(i,2))^2)./v0;
        end
        
        forward.fa.t_shift=t_0-t_ex;
        
    end
    
    % output traveltime
    d{id}=t(:)-forward.fa.t_shift(:);
    
else
    disp(sprintf('%s : forward of type ''%s'' not known',mfilename,forward.type))
    d{id}=[];
    
end




