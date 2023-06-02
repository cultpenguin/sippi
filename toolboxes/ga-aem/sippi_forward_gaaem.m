function [d,forward,prior,data]=sippi_forward_gaaem(m,forward,prior,data,id,im)
% sippi_forward_gaaem: GA-AEM (https://github.com/GeoscienceAustralia/ga-aem) forward response
%
% Call 
%   [d,forward,prior,data]=sippi_forward_gaaem(m,forward,prior,data,id,im)
%
%  forward.thickness [Nm-1,1]
%  forward.isResistivity [0(def),1] is m{1} resisticity or conductivity
%  forward.log10 [0,1(def)] Is m%  forward.log10 [0,1] Is m{1} in log10 space?
%  forward.log10_data [0(def),1] Should output be in log-space?
%
% See: https://github.com/GeoscienceAustralia/ga-aem
%
if nargin==0;
    im=1;
    m{im}=[5.05 0.1 0.05 0.001];
    [d,forward]=sippi_forward_gaaem(m);
    sippi_plot_data_gaaem(d,[],forward);
    return
end


if nargin<2;    forward.null='';end
if nargin<3;    prior{1}.null='';end
if nargin<4;    data{1}.null='';end
if nargin<5;    id=1;end
if nargin<6;
    if isfield(forward,'im')
        im=forward.im;
    else
        im=1;
    end
end

%% 
if size(m{1},2)==1
    m{1}=m{1}';
end



d{1}=1;
%% load library
if ~isfield(forward,'library_loaded')
    forward.library_loaded=0;
end
if forward.library_loaded==0    
    gatdaem1d_loadlibrary()
    forward.library_loaded=1;
end
%% SETUP THE TEM SYSTEM
% get LM
if ~isfield(forward,'LM')
    wdir = fileparts(which('sippi_forward_gaaem.m'));
    forward.LM.stmfile = [wdir,filesep,'stmfiles',filesep,'Skytem-LM.stm'];
end
if ~isfield(forward.LM,'nw')
    forward.LM.hS  = gatdaem1d_getsystemhandle(forward.LM.stmfile);
    forward.LM.nw  = gatdaem1d_nwindows(forward.LM.hS);
    forward.LM.wt  = gatdaem1d_windowtimes(forward.LM.hS);
    forward.LM.wfm = gatdaem1d_waveform(forward.LM.hS);
end
if ~isfield(forward.LM,'hS')
    forward.LM.hS  = gatdaem1d_getsystemhandle(forward.LM.stmfile);
end

% get HM
if ~isfield(forward,'HM')
    wdir = fileparts(which('sippi_forward_gaaem.m'));
    forward.HM.stmfile = [wdir,filesep,'stmfiles',filesep,'Skytem-HM.stm'];
end

if ~isfield(forward.HM,'nw')
    forward.HM.hS  = gatdaem1d_getsystemhandle(forward.HM.stmfile);
    forward.HM.nw  = gatdaem1d_nwindows(forward.HM.hS);
    forward.HM.wt  = gatdaem1d_windowtimes(forward.HM.hS);
    forward.HM.wfm = gatdaem1d_waveform(forward.HM.hS);
end
if ~isfield(forward.HM,'hS')
    forward.HM.hS  = gatdaem1d_getsystemhandle(forward.HM.stmfile);
end
    
if ~isfield(forward,'dtype')
    forward.dtype = gatdaem1d_derivativestruct();
end
if ~isfield(forward,'dlayer')
    forward.dlayer=1;
end

% GEOMETRY
if ~isfield(forward,'G')
    forward.G.null='';
end
if ~isfield(forward.G,'tx_height');forward.G.tx_height = 30;end
if ~isfield(forward.G,'tx_roll');  forward.G.tx_roll   = 0;end
if ~isfield(forward.G,'tx_pitch'); forward.G.tx_pitch  = 0; end
if ~isfield(forward.G,'tx_yaw');   forward.G.tx_yaw    = 0;end
if ~isfield(forward.G,'txrx_dx');  forward.G.txrx_dx   = -12.62;  end
if ~isfield(forward.G,'txrx_dy');  forward.G.txrx_dy   = 0; end
if ~isfield(forward.G,'txrx_dz');  forward.G.txrx_dz   = +2.16;end
if ~isfield(forward.G,'rx_roll');  forward.G.rx_roll   = 0;  end
if ~isfield(forward.G,'rx_pitch'); forward.G.rx_pitch  = 0; end
if ~isfield(forward.G,'rx_yaw');   forward.G.rx_yaw    = 0;end

%% CHECK IF tx_height is set as an input
if ~isfield(forward,'use_tx_height');
    forward.use_tx_height  = 0;
    if nargin>2
        for ip=1:length(prior);
            if strcmp(lower(prior{ip}.name),'tx_height')
                forward.use_tx_height=ip;
            end
        end
    end
end
% update tx_height from m{2}
if forward.use_tx_height>0
    if nargin>2
        %disp(forward.G.tx_height);
        %disp(forward.use_tx_height)
        forward.G.tx_height = m{forward.use_tx_height};
        %disp(forward.G.tx_height);
    end
end

%% SETUP THE MODEL
if ~isfield(forward,'thickness')
    if nargin<3
        forward.thickness = ones(1,length(m{1})-1).*20;
    else
        if length(prior{im}.z)>1
            forward.thickness = diff(prior{im}.z);
        else
            forward.thickness = diff(prior{im}.x);
        end
    end
end
if ~isfield(forward,'log10')
    forward.log10=1;
end
if ~isfield(forward,'log10_data')
    forward.log10_data=0;
end
if ~isfield(forward,'isResistivity')
    forward.isResistivity=0;
end

if (forward.isResistivity)
    m{im}=1./m{im};
end

if (forward.log10==1)
    E.conductivity = 10.^(m{im}');
else
    E.conductivity = m{im}';
end
E.thickness = forward.thickness;

if ~isfield(forward,'compress'); forward.compress=0;end
if forward.compress == 1;
    conductivity = E.conductivity;
    thickness = E.thickness;
    ct = cumsum(thickness);
    
    
   
    ii = find(abs(diff(conductivity))>0);
    if length(ii)==0
        t=[];
    else
        i0=1;
        for i=1:length(ii);
            t(i)=sum(thickness(i0:ii(i)));
            i0=ii(i)+1;
        end
    end
    try
        E.thickness=t;
    catch
        keyboard
    end
    
    E.conductivity=conductivity([1;1+ii(:)]);
    
    
    %cumsum(E.thickness)
    %E.conductivity

end


%% Compute forward response
%LD = gatdaem1d_derivative(forward.LM.hS,forward.dtype.dC,forward.dlayer);
%HD = gatdaem1d_derivative(forward.HM.hS,dtype.dC,forward.dlayer);

useFast=0;
if useFast
    LR = gatdaem1d_forwardmodel(forward.LM.hS,forward.G,E);
    HR = gatdaem1d_forwardmodel(forward.HM.hS,forward.G,E);
    d{id}=-1*[LR.SZ(:);HR.SZ(:)];
else        
    L = gatdaem1d_fm_dlogc(forward.LM.hS,forward.G,E);
    H = gatdaem1d_fm_dlogc(forward.HM.hS,forward.G,E);
    d{id}=-1*[L.FM.SZ(:);H.FM.SZ(:)];
end


if (forward.log10_data==1)
    d{id}=log10(d{id});
end

%if forward.use_tx_height>0
%    d{id+1}=m{2};
%end

% free memory?
if ~isfield(forward,'free_memory');
    forward.free_memory=0;
end
if forward.free_memory==1;
    gatdaem1d_freesystemhandle(forward.LM.hS);
    gatdaem1d_freesystemhandle(forward.HM.hS);
    gatdaem1d_unloadlibrary();
    forward.library_loaded=0;
end

% TMH: 2019-11-14: Always output all data, and use data{1}.i_use, only in
% sippi_likelihood!!!
% if nargin>3
%     if isfield(data{1},'i_use');
%         d{id}=d{id}(data{id}.i_use);
%     end
% end



