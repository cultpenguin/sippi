% sippi_least_squares
%
% Call : 
%    [m_reals,m_est,Cm_est]=sippi_least_squares(data,prior,forward,n_reals,lsq_type,id,im);
%
%
%
%   lsq_type : 'lsq' (def), classical least squares 
%              'error_sim', simulation through error simulation
%              'visim', simulation through SGSIM of DSSIM
%
function [m_reals,m_est,Cm_est,options,forward]=sippi_least_squares(data,prior,forward,n_reals,lsq_type,id,im,options);

if nargin<4, n_reals=15;end
if nargin<5, lsq_type='lsq';end
if nargin<6, id=1;end
if nargin<7, im=1;end
options.null='';
%% CHOOSE NAME
if ~isfield(options,'txt')
    options.txt='sippi_least_squares';
end
options.txt=sprintf('%s_least_squares_%s',datestr(now,'YYYYmmdd_HHMM'),options.txt);

options.txt=sprintf('%s_%s',options.txt,lsq_type);
try
    if forward.linear==1;
        t='linear';
    else
        t='nonlin';
    end
    options.txt=sprintf('%s_%s',options.txt,t);
end

disp(options.txt)


%% MODEL COVARINCE
if ~isfield(prior{im},'Cm');
    prior{im}.Cm=precal_cov([prior{im}.xx(:) prior{im}.yy(:) prior{im}.zz(:)],[prior{im}.xx(:) prior{im}.yy(:) prior{im}.zz(:)],prior{im}.Va);
else
    disp(sprintf('%s : Model covariance could not be set in prior{%d}.Cm',mfilename,im))
end

%% MODEL COVARINCE
if ~isfield(data{id},'CD');
    m=sippi_prior(prior);
    if isfield(forward,'forward_function');
        [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior,data,id,im);
    else
        [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id,im);
    end
   [logL,L,data]=sippi_likelihood(d,data,id);
    if ~isfield(data{id},'CD');
        disp(sprintf('%s : Data covariance could not be found/set in data{%d}.CD',mfilename,id))
    end
end

%% MODEL COVARINCE
if ~isfield(forward,'G');
    try
        m=sippi_prior(prior);
        if isfield(forward,'forward_function');
            [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior,data,id,im);
        else
            [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id,im);
        end
    end
    if ~isfield(forward,'G');
        disp(sprintf('%s : No forward operator G found in forward',mfilename))
    end
end

if ~isfield(data{id},'d0');data{id}.d0=0;end
if ~isfield(prior{im},'m0');prior{im}.m0=0;end


%% READY TO INVERT NOW SETUP INVERSION

i_use=data{1}.i_use;
if length(data{id}.d0)==1
    dd_obs=data{id}.d_obs(i_use)-data{id}.d0;
else
    dd_obs=data{id}.d_obs(i_use)-data{id}.d0(i_use);
end    
if (strcmp(lsq_type,'least_squares')|strcmp(lsq_type,'lsq'));
    % CLASSICAL LEAST SQUARES
    disp(sprintf('%s : solving lsq using ''%s'' type inversion',mfilename,lsq_type))
    [m_est,Cm_est]=least_squares_inversion(forward.G,prior{im}.Cm,data{id}.CD(i_use,i_use),prior{im}.m0,dd_obs);
    %[m_est,Cm_est]=least_squares_inversion(forward.G,prior{im}.Cm,data{id}.CD,prior{im}.m0,data{id}.d_obs-data{id}.d0);
    m_reals=gaussian_simulation_cholesky(m_est,Cm_est,15);
    %[reals,etype_mean,etype_var]=sippi_get_sample(data,prior,id,im,n_reals);

elseif (strcmp(lsq_type,'visim')); 
    disp(sprintf('%s : solving lsq using ''%s'' type inversion',mfilename,lsq_type))
    %V=visim_init(prior{im}.x,prior{im}.y,prior{im}.z)
    x=prior{im}.x;y=prior{im}.y,z=prior{im}.z;
    G=forward.G;
    Cd=data{1}.CD;
    d_obs=data{1}.d_obs;
    m0=prior{im}.m0;
    V=G_to_visim(x,y,z,d_obs,G,Cd,m0);
    V.nsim=n_reals;
    V=visim_set_variogram(prior{im}.Va,V);
    V=visim(V);
    % export data
    m_est=V.etype.mean';
    m_var=V.etype.var';
    
    m_est=m_est(:);
    Cm_est=m_var(:);
    %Cm_est=diag(m_var(:));
    m_reals=zeros(prod(size(V.etype.mean)),n_reals);
    for i=1:n_reals
        if V.nz==1
            m=V.D(:,:,i)';
        else
            m=V.D(:,:,:,i)';
        end
        m_reals(:,i)=m(:);
    end
    
    
    
elseif (strcmp(lsq_type,'error_sim'));
    disp(sprintf('%s : solving lsq using ''%s'' type inversion',mfilename,lsq_type))
    %V=visim_init(prior{im}.x,prior{im}.y,prior{im}.z)
    x=prior{im}.x;y=prior{im}.y,z=prior{im}.z;
    G=forward.G;
    Cd=data{1}.CD;
    d_obs=data{1}.d_obs;
    m0=prior{im}.m0;
    V=G_to_visim(x,y,z,d_obs,G,Cd,m0);
    V.nsim=n_reals;
    V=visim_set_variogram(prior{im}.Va,V);
    
    V=visim_error_sim(V);
    % export data
    m_est=V.etype.mean';
    m_var=V.etype.var';
    
    m_est=m_est(:);
    Cm_est=m_var(:);
    %Cm_est=diag(m_var(:));
    m_reals=zeros(prod(size(V.etype.mean)),n_reals);
    for i=1:n_reals
        if V.nz==1
            m=V.D(:,:,i)';
        else
            m=V.D(:,:,:,i)';
        end
        m_reals(:,i)=m(:);
    end
    
    
else
    disp(sprintf('%s : ''%s'' type inversion not supported',mfilename,lsq_type))
end

%% SCALE M_EST
x=prior{im}.x;y=prior{im}.y;z=prior{im}.z;
if prior{im}.dim(3)>1
    % 3D
    m_est=reshape(m_est,length(y),length(x),length(z));
elseif prior{im}.dim(2)>1
    % 2D
    m_est=reshape(m_est,length(y),length(x));
else
    % 1D
end



%% EXPORT REALIZATIONS TO DISK
try;
    mkdir(options.txt);
end
% REALS
filename_asc{im}=sprintf('%s%s%s_m%d%s',options.txt,filesep,options.txt,im,'.asc');
fid=fopen(filename_asc{im},'w');
for i=1:n_reals
    fprintf(fid,' %10.7g ',m_reals(:,i));
    fprintf(fid,'\n');
end
fclose(fid);

filename_m_est{im}=sprintf('%s%s%s_m%d_mest%s',options.txt,filesep,options.txt,im,'.asc');
fid=fopen(filename_m_est{im},'w');
fprintf(fid,' %10.7g ',m_est(:));
fclose(fid);

filename_Cm_est{im}=sprintf('%s%s%s_m%d_Cmest%s',options.txt,filesep,options.txt,im,'.asc');
fid=fopen(filename_Cm_est{im},'w');
fprintf(fid,' %10.7g ',Cm_est(:));
fclose(fid);

filename_mat=sprintf('%s%s%s.mat',options.txt,filesep,options.txt);
save(filename_mat);
