function [reals_mat,etype_mean,etype_var,reals_all,ite_num]=sippi_get_sample(im,n_reals,skip_seq_gibbs,data,prior,options);
% sippi_get_sample: Get a posterior sample
%
% Call :
%  [reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(im,n_reals,skip_seq_gibbs,data,prior,options);
%
%    im: A priori model type
%    n_reals: Number of realizations to return
%    skip_seq_gibbs [1] Skip all realization where sequential gibbs is enabled
%                   [0] Use all realization
%    data: SIPPI data structure
%    prior: SIPPI prior structure
%    options: options structure when running sippi_metropolis
%
%
% If located in a SIPPI output folder one can simple use :
%    [reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(im,n_reals);
% or
%    skip_seq_gibbs=0;
%    [reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(im,n_reals,skip_seq_gibbs);
%
%


if nargin<4;
    im_chosen=im;
    % LOAD FROM MAT FILES
    [p,matfile]=fileparts(pwd);
    load(matfile);
    im=im_chosen;
end


if ~exist('n_reals','var');
    n_reals=15;
end
if ~exist('id','var');id=1;end
if ~exist('im','var');im=1;end
if ~exist('options','var');
    options.null='';
end
if ~isfield(options,'mcmc');
    options.mcmc.null='';
end
if ~isfield(options.mcmc,'i_sample');
    skip_seq_gibbs=0;
    options.mcmc.i_sample=1;
end


if ~exist('skip_seq_gibbs','var');
    skip_seq_gibbs=1; % only consider posterior realization AFTER seq gibbs has finished
    %skip_seq_gibbs=0; % coniser posterior realization from iteration number 1
end

x=prior{im}.x;y=prior{im}.y;z=prior{im}.z;


%% BUG/19062014 : m_est should go in forward structure
if exist('m_est','var')|isfield(options,'m_est');
    if isfield(options,'m_est');
        m_est=options.m_est;
        Cm_est=options.Cm_est;
    end
        
    % LEAST SQUARES TYPE INVERSION
    reals=gaussian_simulation_cholesky(m_est,Cm_est,n_reals)';
    reals_all=reals; % dummy output   
    ite_num=1:1:n_reals;
    etype_mean=m_est;
    etype_var=diag(Cm_est);

    if prior{im}.dim(3)>1
        etype_var=reshape(etype_var,length(y),length(x),length(z));
    elseif prior{im}.dim(2)>1
        etype_var=reshape(etype_var,length(y),length(x));
    end


else
    
    if ~isfield(options,'txt');
        [p,options.txt]=fileparts(pwd);
    end
    try
        reals=load(sprintf('%s_m%d.asc',options.txt,im));
    catch
        reals=load(sprintf('%s%s%s_m%d.asc',options.txt,filesep,options.txt,im));
    end
    
    n_reals=min([n_reals,size(reals,1)]);
    
    reals_all=reals;
    n_reals_all=size(reals_all,1);
    
    % GET ITERATION NUMBER AFTER SEQ GIBBS HAS FINISHED, i1
    try
        i1_post=prior{im}.seq_gibbs.i_update_step_max./options.mcmc.i_sample;
    catch
        i1_post=1;
    end
    i1_post=max([i1_post 1]);
    
    if skip_seq_gibbs==0;
        i1_post=1;
    end
    
    %% TAKE OUT ALL 'POSTERIOR' REALIZATION AND THE ITERATION NUMBER
    nr=size(reals,1);
    
    % all posterior sampels
    ii_post_reals=ceil(i1_post./options.mcmc.i_sample):1:n_reals_all;
    ni_post_reals=ii_post_reals*options.mcmc.i_sample;
    
    reals=reals(ii_post_reals,:);
    nr=size(reals,1);
    
    %% GET ETYPES
    etype_mean=mean(reals);
    etype_var=var(reals);
    if prior{im}.dim(3)>1
        etype_mean=reshape(etype_mean,length(y),length(x),length(z));
        etype_var=reshape(etype_var,length(y),length(x),length(z));
    elseif prior{im}.dim(2)>1
        etype_mean=reshape(etype_mean,length(y),length(x));
        etype_var=reshape(etype_var,length(y),length(x));
    end
    
    %% TAKE OUT ONLY n_reals REALIZATION
    if (nr<n_reals)
        n_reals=nr;
    end
    N=length(ii_post_reals);
    
    i_use = ceil(linspace(1,N,n_reals));
    
    ii=ii_post_reals(i_use);
    ni=ni_post_reals(i_use);
    
    
    n_reals=min([n_reals size(reals,1)]);
    
    ii=ceil(linspace(1,size(reals,1),n_reals));
    ite_num=ii*options.mcmc.i_sample;
    reals=reals(i_use,:);
end

%% GET REQUESTED SAMPLE
if n_reals<1
    reals_mat=[];
    disp(sprintf('%s : Number of ''realizations'' is less than one!',mfilename));
end

for i=1:n_reals
    if prior{im}.dim(3)>1
        % 3D
        reals_mat(:,:,:,i)=reshape(reals(i,:),length(y),length(x),length(z));
    elseif prior{im}.dim(2)>1
        % 2D
        reals_mat(:,:,i)=reshape(reals(i,:),length(y),length(x));
    else
        % 1D
        reals_mat(:,i)=reals(i,:);
    end
end

