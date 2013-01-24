function [reals_mat,etype_mean,etype_var,reals]=sippi_get_sample(data,prior,id,im,n_reals,options);
% sippi_get_sample Get a posterior sample 
%
% Call :
%  [reals,etype_mean,etype_var]=sippi_get_sample(data,prior,id,im,n_reals,options);
%
if ~exist('n_reals','var');
    n_reals=15;
end
if ~exist('id','var');id=1;end
if ~exist('im','var');im=1;end
if ~exist('options','var');options.null='';end
%% LSQ
if isfield(data{id},'m_est');
    reals=gaussian_simulation_cholesky(data{im}.m_est,data{im}.Cm_est,n_reals)';
    etype_mean=data{id}.m_est;
    etype_var=diag(data{id}.Cm_est);
else
    
    if ~isfield(options,'txt');
        [p,options.txt]=fileparts(pwd);
    end
    try        
        reals=load(sprintf('%s_m%d.asc',options.txt,im));
    catch
        reals=load(sprintf('%s%s%s_m%d.asc',options.txt,filesep,options.txt,im));
    end
    try
        i1_post=prior{im}.seq_gibbs.i_update_step_max./options.mcmc.i_sample;
    catch
        i1_post=1;
    end
    
    i1_post=max([i1_post 1]);
    
    nr=size(reals,1);
    if (nr<n_reals)
        n_reals=nr;
    end
    i1=ceil(size(reals,1)/n_reals);
    i1=max([i1 i1_post]);
    
    % ONLY REALS FROM POSTERIOR
    if i1_post>n_reals
        i1_post=1;
    end
    reals=reals(i1_post:end,:);
    
    ii=ceil(linspace(1,size(reals,1),n_reals));
    
    
    x=prior{im}.x;y=prior{im}.y;z=prior{im}.z;
   
  
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
    
    
    %% GET REQUESTED SAMPLE
    
        
        for i=1:n_reals
            if prior{im}.dim(3)>1
                % 3D
                reals_mat(:,:,:,i)=reshape(reals(ii(i),:),length(y),length(x),length(z));
            elseif prior{im}.dim(2)>1
                % 2D
                reals_mat(:,:,i)=reshape(reals(ii(i),:),length(y),length(x));
            else
                % 1D
                reals_mat(:,i)=reals(ii(i),:);
            end
        end
        
end
