% sippi_least_squares Least squares type inversion for SIPPI
%
% Call :
%   [m_est,Cm_est,m_reals,options,data,prior,forward]=sippi_least_squares(data,prior,forward,options);
%
%   options.lsq.type    : LSQ type to use ('lsq' (classical linear leqast squares) is the default)
%   options.lsq.n_reals : Number of realizations to generate
%   options.lsq.plot : [0/1] show figures or not def->0.
%   options.lsq.save_data : [0/1] save realizations to output folder. def->1.
%
%
% TMH/01/2017
%
% See also sippi_rejection, sippi_metropolis
%



%              'error_sim', simulation through error simulation
%              'visim', simulation through SGSIM of DSSIM
%
function [m_est,Cm_est,m_reals,options,data,prior,forward]=sippi_least_squares(data,prior,forward,options);
id=1;
im=1;

m_reals=[];
m_est=[];
Cm_est=[];
options.lsq.null=''; % make sutre options.lsq exists

% number of realizations
if ~isfield(options.lsq,'n_reals');
    options.lsq.n_reals=50;
end
if ~isfield(options.lsq,'type');
    options.lsq.type='lsq';
    %options.lsq.type='error_sim';
    %options.lsq.type='visim';
end

% compute reals
if ~isfield(options.lsq,'compute_reals')
    if nargout>2
        options.lsq.compute_reals=1;
    else
        options.lsq.compute_reals=0;
        options.lsq.save_data=0;
    end
end

% save data?
if ~isfield(options.lsq,'save_data')
    options.lsq.save_data=options.lsq.compute_reals;
end


% plot?
if ~isfield(options.lsq,'plot')
    options.lsq.plot=0;
end


%% CHOOSE NAME
if ~isfield(options,'txt')
    options.txt=mfilename;%'sippi_least_squares';
end


options.txt=sprintf('%s_%s_%s',datestr(now,'YYYYmmdd_HHMM'),options.txt,options.lsq.type);
sippi_verbose(sprintf('%s: output folder: %s ',mfilename,options.txt),1)

%% MODEL COVARINCE
if ~isfield(options.lsq,'Cm');
    prior=sippi_prior_init(prior);
    if isfield(prior{im},'Cmat');
        options.lsq.Cm=prior{im}.Cmat;
    else
        prior=sippi_prior_init(prior);
        options.lsq.Cm=precal_cov([prior{im}.xx(:) prior{im}.yy(:) prior{im}.zz(:)],[prior{im}.xx(:) prior{im}.yy(:) prior{im}.zz(:)],prior{im}.Va);
    end
end
if ~isfield(options.lsq,'Cm');
    sippi_verbose(sprintf('%s: Could not model covariance Cm. Please use a Gaussian prior or set prior{%d}.Cmat',mfilename,id),0)
else
    sippi_verbose(sprintf('%s: Model covariance, Cm, set in options.lsq.Cm',mfilename),1)
end

%% DATA COVARINCE
if isfield(data{id},'CD');
    options.lsq.Cd=data{id}.CD;
else
    % solve the forwrad problem and make use of data{1}.CD if it exists
    m=sippi_prior(prior);
    [d,forward,prior,data]=sippi_forward(m,forward,prior,data);
    [logL,L,data]=sippi_likelihood(d,data,id);
    try
        options.lsq.Cd=data{id}.CD;
    end
    try
        % BIAS?
        options.lsq.d0=data{id}.dt;
    end
end

if ~isfield(options.lsq,'Cd');

    if isfield(data{id},'d_std');
        if length(data{id}.d_std)==1;
            options.lsq.Cd=eye(length(data{id}.d_obs)).*data{id}.d_std.^2;
        else
            options.lsq.Cd=diag(data{1}.d_std.^2);
        end
    end
    if isfield(data{id},'d_var');
        if length(data{id}.d_var)==1;
            options.lsq.Cd=eye(length(data{id}.d_obs)).*data{id}.d_var;
        else
            options.lsq.Cd=diag(data{1}.d_var);
        end
    end
end

if ~isfield(options.lsq,'Cd');
    sippi_verbose(sprintf('%s: Could not data covariance Cd. Please use a Gaussian noise model in data{%d}',mfilename,id),0)
else
    sippi_verbose(sprintf('%s: Data covariance, Cd, set in options.lsq.Cd',mfilename),1)
end


%% CHECK FOR FORWARD OPERATOR
if ~isfield(forward,'G');
    % assume the forward operator is output in forward.G if sippi_forward
    % is run
    try
        m=sippi_prior(prior);
        if isfield(forward,'forward_function');
            [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior,data,id,im);
        else
            [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id,im);
        end
    end
    if ~isfield(forward,'G');
        sippi_verbose(sprintf('%s : No forward operator G found in forward',mfilename),0)
    end
end
try
    options.lsq.G=forward.G;
end
if ~isfield(options.lsq,'G');
    sippi_verbose(sprintf('%s: linear forward operator G is not set. Please set in in forward.G',mfilename,id),0)
else
    sippi_verbose(sprintf('%s: linear forward operator G set in options.lsq.G',mfilename),1)
end

%% M
if ~isfield(prior{im},'m0');
    prior{im}.m0=0;    
end

if length(prior{im}.m0)==1;
    nm=size(options.lsq.Cm,1);
    options.lsq.m0=ones(nm,1).*prior{im}.m0;
else
    options.lsq.m0=prior{im}.m0(:);
end


if isfield(forward,'linear_m');
    if length(forward.linear_m)==1;
        nm=size(options.lsq.Cm,1);
        options.lsq.m0=ones(nm,1).*forward.linear_m;
    else
        options.lsq.m0=forward.linear_m;
    end
end

sippi_verbose(sprintf('%s: setting options.lsq.m0=prior{%d}.m0',mfilename,im),1)

%% D
if ~isfield(options.lsq,'d0');
    options.lsq.d0=data{id}.d_obs.*0;
    sippi_verbose(sprintf('%s: setting options.lsq.d0=0;',mfilename),1)
end
options.lsq.d_obs=data{id}.d_obs;
sippi_verbose(sprintf('%s: setting options.lsq.d_obs=data{%d}.d_obs',mfilename,id),1)

%%
i_use=data{1}.i_use;
n_use=length(i_use);
sippi_verbose(sprintf('%s : Linear least squares using ''%s'' type inversion',mfilename,options.lsq.type),1)
if (strcmp(options.lsq.type,'least_squares')|strcmp(options.lsq.type,'lsq'));
    % CLASSICAL LEAST SQUARES
    if  n_use==size(options.lsq.G,1);
        [m_est{im},Cm_est{im}]=least_squares_inversion(options.lsq.G,options.lsq.Cm,options.lsq.Cd(i_use,i_use),options.lsq.m0,options.lsq.d_obs(i_use)-options.lsq.d0(i_use));
    else
        [m_est{im},Cm_est{im}]=least_squares_inversion(options.lsq.G(i_use,:),options.lsq.Cm,options.lsq.Cd(i_use,i_use),options.lsq.m0,options.lsq.d_obs(i_use)-options.lsq.d0(i_use));
    end
    if (options.lsq.compute_reals==1)
        m_reals=gaussian_simulation_cholesky(m_est{im},Cm_est{im},options.lsq.n_reals);
    end
    % elseif (strcmp(lsq_type,'visim'));
    %     %% LSQ USING SEQUENTIAL SIMULATION IN VISIM
    %     sippi_verbose(sprintf('%s : solving lsq using ''%s'' type inversion',mfilename,lsq_type))
    %     %V=visim_init(prior{im}.x,prior{im}.y,prior{im}.z)
    %     x=prior{im}.x;y=prior{im}.y,z=prior{im}.z;
    %     G=forward.G;
    %     Cd=data{1}.CD;
    %     d_obs=data{1}.d_obs;
    %     m0=prior{im}.m0;
    %     V=G_to_visim(x,y,z,d_obs,G,Cd,m0);
    %     V.nsim=n_reals;
    %     V=visim_set_variogram(prior{im}.Va,V);
    %     V=visim(V);
    %     % export data
    %     m_est=V.etype.mean';
    %     m_var=V.etype.var';
    %
    %     m_est=m_est(:);
    %     Cm_est=m_var(:);
    %     %Cm_est=diag(m_var(:));
    %     m_reals=zeros(prod(size(V.etype.mean)),n_reals);
    %     for i=1:n_reals
    %         if V.nz==1
    %             m=V.D(:,:,i)';
    %         else
    %             m=V.D(:,:,:,i)';
    %         end
    %         m_reals(:,i)=m(:);
    %     end
    % elseif (strcmp(lsq_type,'error_sim'));
    %     %% LSQ USING ERRROR SIMULATION IN VISIM
    %     sippi_verbose(sprintf('%s : solving lsq using ''%s'' type inversion',mfilename,lsq_type))
    %     %V=visim_init(prior{im}.x,prior{im}.y,prior{im}.z)
    %     x=prior{im}.x;y=prior{im}.y,z=prior{im}.z;
    %     G=forward.G;
    %     Cd=data{1}.CD;
    %     d_obs=data{1}.d_obs;
    %     m0=prior{im}.m0;
    %     V=G_to_visim(x,y,z,d_obs,G,Cd,m0);
    %     V.nsim=n_reals;
    %     V=visim_set_variogram(prior{im}.Va,V);
    %
    %     V=visim_error_sim(V);
    %     % export data
    %     m_est=V.etype.mean';
    %     m_var=V.etype.var';
    %
    %     m_est=m_est(:);
    %     Cm_est=m_var(:);
    %     %Cm_est=diag(m_var(:));
    %     m_reals=zeros(prod(size(V.etype.mean)),n_reals);
    %     for i=1:n_reals
    %         if V.nz==1
    %             m=V.D(:,:,i)';
    %         else
    %             m=V.D(:,:,:,i)';
    %         end
    %         m_reals(:,i)=m(:);
    %     end
    %
    %
else
    sippi_verbose(sprintf('%s : ''%s'' type inversion not supported',mfilename,lsq_type),-1)
end

%% SCALE M_EST
x=prior{im}.x;y=prior{im}.y;z=prior{im}.z;
if prior{im}.dim(3)>1
    % 3D
    m_est{im}=reshape(m_est{im},length(y),length(x),length(z));

elseif prior{im}.dim(2)>1
    % 2D
    m_est{im}=reshape(m_est{im},length(y),length(x));
    Cm_est_diag = reshape(diag(Cm_est{im}),length(y),length(x));
else
    % 1D
    Cm_est_diag = diag(Cm_est{im});
end


%% EXPORT REALIZATIONS TO DISK
if (options.lsq.save_data==1)

    options.cwd=pwd;
    try;
        mkdir(options.txt);
    end
    % REALS
    filename_asc{im}=sprintf('%s%s%s_m%d%s',options.txt,filesep,options.txt,im,'.asc');
    fid=fopen(filename_asc{im},'w');
    for i=1:options.lsq.n_reals
        fprintf(fid,' %10.7g ',m_reals(:,i));
        fprintf(fid,'\n');
    end
    fclose(fid);


    filename_m_est{im}=sprintf('%s%s%s_m%d_mest%s',options.txt,filesep,options.txt,im,'.asc');
    sippi_verbose(sprintf('%s: Writing m_est to %s',mfilename,filename_m_est{im}),1);
    %fid=fopen(filename_m_est{im},'w');fprintf(fid,' %10.7g ',m_est(:));fclose(fid);

    filename_Cm_est{im}=sprintf('%s%s%s_m%d_Cmest%s',options.txt,filesep,options.txt,im,'.asc');
    sippi_verbose(sprintf('%s: Writing Cm_est to %s',mfilename,filename_Cm_est{im}),1);
    %fid=fopen(filename_Cm_est{im},'w');fprintf(fid,' %10.7g ',Cm_est(:));fclose(fid);

    filename_mat=sprintf('%s%s%s.mat',options.txt,filesep,options.txt);
    sippi_verbose(sprintf('%s: Writing %s',mfilename,filename_mat),1);
    save(filename_mat);
end

%% PLOT
if options.lsq.plot==1;
    figure(71);clf;
    subplot(1,2,1);
    m{1}=m_est;
    sippi_plot_prior(prior,m,im,0,gca);
    colorbar
    subplot(1,2,2);
    m{1}=Cm_est_diag;
    sippi_plot_prior(prior,m,im,0,gca);
    if isfield(prior{im},'cax_var');
        caxis(prior{im}.cax_var);
    else
        caxis([0 max(m{1}(:))])
    end
    colorbar
    filename_png=sprintf('%s%s%s_mEst_CmEst.mat',options.txt,filesep,options.txt);
    print_mul(filename_png);
end


