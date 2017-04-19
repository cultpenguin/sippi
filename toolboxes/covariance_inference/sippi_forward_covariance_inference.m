% sippi_forward_covariance_inference : Probabilitsic covariance inference
%
% Call :
%  [d,forward,prior,data]=sippi_forward_covariance_inference(m,forward,prior,data,id,im)
%
%
%  forward.pos_known : [x' y' z'],  [ndata,ndim] with position of observed data
%  forward.G : Forward operator
%
%  Prior covariance model, N(m0,Cm) is chosen as
%  forward.m0 : initial mean model
%  forward.Cm/forward.Va : inital covariance model
%
%  prior{im}.m0
%  prior{im}.type (round(type) is itype)
%  prior{im}.range_1
%  prior{im}.range_2
%  prior{im}.range_3
%  prior{im}.ang_1
%  prior{im}.ang_2
%  prior{im}.ang_3
%  prior{im}.sill
%  prior{im}.nugget_fraction
%
% See Hansen et al. (2014). A general probabilistic approach for inference of Gaussian model parameters from noisy data of point and volume support. 
% Mathematical Geosciences 47(7), pp 843-865. published online 09-2014. doi:10.1007/s11004-014-9567-5 
%
% See also: sippi_forward
%
function [d,forward,prior,data]=sippi_forward_covariance_inference(m,forward,prior,data,id,im)

if nargin<6,
    try
        im=length(prior);
    catch
        im=1;
    end
end
if nargin<5,
    try
        id=length(data);
    catch
        id=1;
    end
end

if nargin<4, data{1}.null='';end

if ~isfield(forward,'Va');
    if isfield(forward,'Cm');
        forward.Va=forward.Cm;
    else
        disp(sprintf('Initial covariance model MUST be set',mfilename));
    end
end
if (isstr(forward.Va))
    forward.Va=deformat_variogram(forward.Va);
end

ns=length(forward.Va);

if ~isfield(forward,'stabilize');
    forward.stabilize=0;
end

if ~isfield(forward,'point_support');
    forward.point_support=0;
end

if ~isfield(data{1},'i_use')
    data{1}.i_use=1:1:length(data{1}.d_obs);
end


% make sure the nugget is present.
if ns==1;
    forward.Va(2)=forward.Va(1);
    forward.Va(1).type='Nug';
    forward.Va(1).par1=0;
    forward.Va(1).par2=0;
    forward.Va(1).itype=0;
    ns=2;
end


%% SET THE COVARIANCE MODEL ACCORDING TO THE CURRENT MODEL

% check dimension
dim=1;
if length(forward.Va(2).par2)==1;dim=1;end
if length(forward.Va(2).par2)==3;dim=2;end
if length(forward.Va(2).par2)==6;dim=3;end

for im=1:length(prior)
    if strcmp(prior{im}.name,'range_1');dim=max([dim 1]);end
    if strcmp(prior{im}.name,'range_2');dim=max([dim 2]);end
    if strcmp(prior{im}.name,'range_3');dim=max([dim 3]);end
    if strcmp(prior{im}.name,'ang_1');dim=max([dim 2]);end
    if strcmp(prior{im}.name,'ang_2');dim=max([dim 3]);end
    if strcmp(prior{im}.name,'ang_3');dim=max([dim 3]);end
end

if dim>1;
    par2=forward.Va(2).par2;
    n_par2=length(par2);
    
    if dim==2
        if n_par2<3
            forward.Va(2).par2=zeros(1,3);
            forward.Va(2).par2(3)=1;
            forward.Va(2).par2(1:n_par2)=n_par2;
        end
    else
        forward.Va(2).par2=zeros(1,6);
        forward.Va(2).par2(5)=1;
        forward.Va(2).par2(6)=1;
        if n_par2==1;
            forward.Va(2).par2(1)=par2(1);
        else
            forward.Va(2).par2(1:2)=par2(1:2);
            forward.Va(2).par2(5)=par2(3);
        end
    end
end

% type
for im=1:length(prior)
    if strcmp(prior{im}.name,'type');
        itype=Speh(m{im});
        itype=max([itype 1]);
        itype=min([itype 3]);
        m{im}=itype;
        forward.Va(ns).itype=itype;
        if forward.Va(ns).itype==1;
            forward.Va(ns).type='Sph';
        elseif forward.Va(ns).itype==2;
            forward.Va(ns).type='Exp';
        elseif forward.Va(ns).itype==3;
            forward.Va(ns).type='Gau';
        elseif forward.Va(ns).itype==0;
            forward.Va(ns).type='Nug';
        end
    end
end


% m0
for im=1:length(prior)
    if strcmp(prior{im}.name,'m0');
        forward.m0=m{im};
    end
end

% range 1
for im=1:length(prior)
    if strcmp(prior{im}.name,'range_1');
        forward.Va(ns).par2(1)=m{im};
    end
end

% range 2
for im=1:length(prior)
    if strcmp(prior{im}.name,'range_2');
        aniso = m{im}./forward.Va(ns).par2(1);
        if dim==2
            forward.Va(ns).par2(3)=aniso;
        elseif dim==3
            forward.Va(ns).par2(5)=aniso;
        end
        
    end
end

% range 3
for im=1:length(prior)
    if strcmp(prior{im}.name,'range_3');
        aniso = m{im}./forward.Va(ns).par2(1);
        forward.Va(ns).par2(6)=aniso;
    end
end

% angle1
for im=1:length(prior)
    if strcmp(prior{im}.name,'ang_1');
        if length(forward.Va(ns).par2)==1;
            % choose no anisotropy unless specified -> range_1 = range_2
            forward.Va(ns).par2(3)=1;
        end
        forward.Va(ns).par2(2)=m{im};
    end
end

% angle2
for im=1:length(prior)
    if strcmp(prior{im}.name,'ang_2');
        forward.Va(ns).par2(3)=m{im};
    end
end

% angle3
for im=1:length(prior)
    if strcmp(prior{im}.name,'ang_3');
        forward.Va(ns).par2(4)=m{im};
    end
end

% sill fraction
for im=1:length(prior)
    sill=sum([forward.Va.par1]);
    nugget_fraction=forward.Va(1).par1(1)/sill;
    if strcmp(prior{im}.name,'sill');
        forward.Va(1).par1(1) = nugget_fraction*m{im};
        forward.Va(2).par1(1) = (1-nugget_fraction)*m{im};
    end
end

% nugget fraction
for im=1:length(prior)
    sill=sum([forward.Va.par1]);
    if strcmp(prior{im}.name,'nugget_fraction');
        forward.Va(1).par1(1) = sill*m{im};
        forward.Va(2).par1(1) = sill*(1-m{im});
    end
end
%
%
% if dim==3;
%     if forward.Va(2).par2(6)==0;forward.Va(2).par2(6)=1;end
%     if forward.Va(2).par2(5)==0;forward.Va(2).par2(5)=1;end
% end
% if dim==2;
%     if forward.Va(2).par2(3)==0;forward.Va(2).par2(3)=1;end
% end

sippi_verbose(sprintf('%s: forward.Va=%s',mfilename,format_variogram(forward.Va)),1)
%%
% make sure CD is recomputed each time
data{1}.recomputeCD=1;
data{1}.full_likelihood=1;

if forward.point_support==1;
    
    if ~isfield(forward,'G');
        forward.G=eye(size(forward.pos_known,1));
    end
    
    % MAKE USE OF FASTER APPROACH FOR POINT SUPPORT
    %if isfield(data{1},'i_use')
    %    i_use=data{1}.i_use;
    %else
        i_use=find(sum(forward.G)); % This line is udefulle when forward.x,.. is set
    %end
    % i.e. when the number of 'model
    % parameters' is larger than the number of
    % 'data
    
    forward.Ct=precal_cov(forward.pos_known(i_use,:),forward.pos_known(i_use,:),forward.Va);
    %forward.Ct=forward.Cm;
    
else
    % VOLUME SUPPORT -> FULL POINT COVARIANCE MUST BE COMPUTED
    
    cal_cm=1;
    if (isfield(forward,'z'));
        % ONLY DO THIS WHEN NMODEL>NDATA (WHEN forward.x,..y,..z is set)
        if length(forward.z)==1;
            % USE FASTER 2D COVARIANCE SETUP
            if ~isfield(forward,'nx'); forward.nx=length(forward.x);end
            if ~isfield(forward,'ny'); forward.ny=length(forward.y);end
            if ~isfield(forward,'dx'); forward.dx=forward.x(2)-forward.x(1);end
            if ~isfield(forward,'dy'); forward.dy=forward.y(2)-forward.y(1);end
            forward.Cm=precal_cov_2d(forward.nx,forward.ny,forward.dx,forward.dy,forward.Va);
            calc_cm=0;
        end
    end
    if cal_cm==1
        forward.Cm=precal_cov_auto([forward.pos_known],forward.Va);
    end
    %forward.Cm=precal_cov([forward.pos_known],[forward.pos_known],forward.Va);
    if ~isfield(forward,'G');
        forward.G=eye(size(forward.Cm));
    end
    
    % COMPUTE EXPECTED DATA TO DATA COVARIANCE !!!
    forward.Ct=forward.G*forward.Cm*forward.G';
    
end

% Stabilize matrix inversion
if isfield(forward,'stabilize');
    forward.Ct=forward.Ct+forward.stabilize.*eye(size(forward.Ct));
end

if isnan(forward.Ct(1))
    keyboard
end

data{1}.Ct=forward.Ct;
if isfield(forward,'m0');
    data{1}.dt=forward.G*(ones(size(forward.G,2),1)*forward.m0(:));
%    data{1}.dt=forward.G(data{1}.i_use,:)*(ones(size(forward.G,2),1)*forward.m0(:));
end
if isfield(forward,'linear_m');
    data{1}.dt=forward.G*forward.linear_m(:);
end
d{id}=0*data{id}.d_obs(data{1}.i_use);


