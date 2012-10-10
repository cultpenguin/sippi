% sippi_forward Simple forward wrapper for SIPPI
%
% Assumes that the actual forward solver has been defined by
% forward.forward_function 
%
% Call:
%   [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id,im)
%
function [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id,im)

%if nargin<4;    forward.null='';end
if nargin<4;    data{1}.null='';end
if nargin<5;    id=1;end
if nargin<6;    im=1;end

if isfield(forward,'forward_function');
    
    if nargin==2;
        [d,forward,prior,data]=feval(forward.forward_function,m,forward);
    elseif nargin==3
        [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior);
    else
        [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior,data,id,im);
    end
else
    disp(sprintf('%s : No forward_function specified in ''forward'' structure',mfilename))
    d=[];
    
end



if isfield(data{id},'sources')
    forward.sources=data{id}.sources;
end
if isfield(data{id},'receivers')
    forward.receivers=data{id}.receivers;
end

S=forward.sources;
R=forward.receivers;

x=prior{im}.x;
y=prior{im}.y;
z=prior{im}.z;

if ~isfield(data{id},'i_use');
    %data{id}.i_use=1:size(data{id}.d_obs,1);
    data{id}.i_use=1:size(forward.sources,1);
end

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
        d{id}=eikonal_traveltime(prior{im}.x,prior{im}.y,prior{im}.z,1./m{im},forward.sources,forward.receivers,data{id}.i_use,forward.eikonal_type);
    else
        % VELOCITY
        d{id}=eikonal_traveltime(prior{im}.x,prior{im}.y,prior{im}.z,m{im},forward.sources,forward.receivers,data{id}.i_use,forward.eikonal_type);
    end
elseif (strcmp(forward.type,'ray')|strcmp(forward.type,'fat'));
    % RAY APPROXIMATION    
    if ~isfield(forward,'linear');forward.linear=0;end
    if (~isfield(forward,'G')|forward.linear==0);
        % ONLY COMPUTE KERNEL IF linear=1 OR IF G DOES NOT EXIST
        if ~isfield(forward,'freq'), forward.freq = 0.1;end
        T=1./ forward.freq;
        
        if ~isfield(forward,'linear_m'),
            forward.linear_m = m{im}.*0+prior{im}.m0;
        end
        if forward.linear==1
            m_use=forward.linear_m;
        else
            m_use=m{im};
        end
        
        if forward.is_slowness==1;
            [K,RAY,Gk,Gray]=tomography_kernel(1./m_use,x,y,z,S(data{id}.i_use,:),R(data{id}.i_use,:),T);
        else
            [K,RAY,Gk,Gray]=tomography_kernel(m_use,x,y,z,S(data{id}.i_use,:),R(data{id}.i_use,:),T);
        end
    
        if strcmp(forward.type,'ray')
            forward.G=Gray;
        else
            forward.G=Gk;
        end
    end
    if forward.is_slowness==1
        d{id}=forward.G*m{im}(:);
    else
        s=1./m{im}(:);
        d{id}=forward.G*s;
    
    end
    
elseif strcmp(forward.type,'born');
    % BORN APPROXIMATION
    % LIU ET AL (2009) for seismic waves
    % Buursink et al () for electromagnetic waves
     if ~isfield(forward,'linear_m'),
            forward.linear_m = m{im}.*0+prior{im}.m0;
     end
     if forward.linear==1
        m_use=forward.linear_m;
    else
        m_use=m{im};
    end
    forward.G=zeros(length(data{id}.i_use),prod(prior{im}.dim));
    j=0;
    for i=data{id}.i_use;%size(forward.sources,1)
        j=j+1;
        progress_txt(i,length(data{id}.d_obs),'setting up kernel');
        [kernel,L,L1_all,L2_all]=kernel_buursink_2d(m_use,x,y,forward.sources(i,:),forward.receivers(i,:),forward.freq);
        forward.G(j,:)=kernel(:);
    end
    if forward.is_slowness==1
        d{id}=forward.G*m{im}(:);
    else
        s=1./m{im}(:);
        d{id}=forward.G*s;
    end
    
elseif strcmp(forward.type,'fw_firstarrival');
    %% ERNST FW CODE WITH KNUDS PARRALLEL MATLAB INTERFACE
    [d,forward,prior,data]=sippi_forward_gpr_ernst_knud(m,forward,prior,data,id,im);
else
    disp(sprintf('%s : forward of type ''%s'' not known',mfilename,forward.type))
    d{id}=[];
    
end
    

%%[logL,L,data]=sippi_likelihood(d,data,id);



