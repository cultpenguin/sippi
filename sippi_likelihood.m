% sippi_likelihood Compute likelihood given an observed dataset
%
% Call
%   [logL,LogL_all,data]=sippi_likelihood(d,data);
%
%
%  data{1}.d_obs [N_data,1] N_data data observations
%  data{1}.d_std [N_data,1] N_data uncorrelated Gaussian STD
%
%  data{1}.d_var [N_data,1] N_data uncorrelated Gaussian variances
%
%
% Gaussian modelization error, N(dt,Ct), is specified as
%  data{1}.dt [N_data,1] : Bias/mean of modelization error
%  data{1}.Ct [N_data,N_data] : Covariance of modelization error
%
%  data{1}.Ct [1,1] : Constant Covariance of modelization error
%                     imples data{1}.Ct=ones(N_data.N_data)*data{1}.Ct;
%
% data{id}.recomputeCD [default=0], if '1' then data{1}.iCD is recomputed
% each time sippi_likelihood is called. This should be used if the noise model
% changes between each call to sippi_likelihood.
%
%  data{id}.full_likelihood [default=]0; if '1' the the full likelihood
%  (including the determinant) is computed. This not needed if the data
%  civariance is constant, but if it changes, then use
%  data{id}.full_likelihood=1;
%
%
%  A new type of noise model can be used as long as it is available in a
%  m file staring with 'sippi_likelihood_'. Further, it should provide the
%  inputs and outputs as sippi_likelihood.m
%  If a noise model has been implemented in the m-files
%  sippi_likelihood_other.m
%  then this can be used to evaluate the likelhood in sippi using
%  data{1}.noise_model='sippi_likelihood_other',
%
%  %% Using multiple data
%  for each data one can specify whether it is used. For example
%  if using three data types (length(data)=3), one can ignore the second
%  type of data (when computing the log-likelihood) using:
%    data{1}.use=1;
%    data{2}.use=0;
%    data{3}.use=1;
%
% % OUTPUT 
%     'logL': sum of all all log-likelihoods
%     'logL_all': Iniduvual log-likelihood for each type of data
%     'data': update data structre 
%
%



%% MAKE SURE Cd CT Ct is robust to simple noise models
%% MAKE sure uncorrelated noise is accounted for in fast way!
%% CHECK THAT GENERALIZED GAUSSIAN WORKS FOR data{1}.d_std. and/or data{1}.d_var !!
%
%
function [logL,L,data]=sippi_likelihood(d,data,id_array)

% Select which data for which likelihood should be computed!!
% This is perhaps a bit slow!!
if  nargin<3
    if ~isfield(data{1},'is_checked');
        for id=1:length(data)
            if ~isfield(data{id},'use');
                data{id}.use=1;
            end
        end
    end
    id_array=[];
    for id=1:length(data)
        if data{id}.use==1;
            id_array=[id_array id];
        end
    end
    % OLD
    % id_array=1; % USE ONLY FIRST DATA
    %id_array=1:length(d); %% USE ALL DATA
    data{1}.is_checked=1;    
end
%id_array=1; % USE ONLY FIRST DATA

logL=zeros(1,length(d));
L=zeros(1,length(d));
for id=id_array;  
    % Check whether 'full_likelihood' is set
    if ~isfield(data{id},'full_likelihood');
        data{id}.full_likelihood=0;
    end
    
    if ~isfield(data{id},'recomputeCD')
        data{id}.recomputeCD=0;
    end
    
    % check whether to use log of data before comptuting misfit?
    if ~isfield(data{id},'use_log')
        data{id}.use_log=0;
    end

%     % Check whether to use user supplied noise model.
%     if isfield(data{id},'noise_model')
%         % next line may be slow...
%         if strfind(data{id}.noise_model,'sippi_likelihood_');
%             [logL,logL_all,data]=feval(data{id}.noise_model,d,data,id);
%             break;
%         end
%     end
%     
    
    if ~isfield(data{id},'noise_model');
        data{id}.noise_model='gaussian';
        %data{id}.noise_model='generalized_gaussian';
        %data{id}.noise_model='laplace';
    end
    
    if ~isfield(data{id},'noise_uncorr');
        %if (~isfield(data{id},'Cd')|~isfield(data{id},'CD'))
        %    % Force uncorrelated noise in case Cd ot CD is not set!!
        %    data{id}.noise_uncorr=1;
        %end
        if (isfield(data{id},'Ct')||isfield(data{id},'Cd')||isfield(data{id},'CD'))
            % Force Correlated noise in case Cd ot CD is not set!!
            data{id}.noise_uncorr=0;
        else
            data{id}.noise_uncorr=1;
        end
        
        
    end
    
    if ~isfield(data{id},'i_use'); data{id}.i_use=1:1:length(data{id}.d_obs);end
    
    
    if strcmp(data{id}.noise_model,'gaussian')&&(data{id}.noise_uncorr==1)
        % UNCORRELATED GAUSSIAN NOISE
        
        % dd=data{id}.d_obs-d{id};
        % d_std could be an array of lenth(data{id}.d_obs)...
        
        if data{id}.use_log==1;
            dd=log(data{id}.d_obs(data{id}.i_use))-log(d{id}(data{id}.i_use));
        else
            dd=data{id}.d_obs(data{id}.i_use)-d{id}(data{id}.i_use);            
        end
        if isfield(data{id},'n_score');
            keyboard
            dd = nscore_mul(dd,data{id}.n_score(data{id}.i_use)); 
        end
        
        if length(data{id}.d_std)==1
            logL(id)=-.5*sum(sum(sum(dd.^2./(data{id}.d_std.^2))));
        else
            logL(id)=-.5*sum(sum(sum(dd.^2./(data{id}.d_std(data{id}.i_use).^2))));
        end
        if data{id}.full_likelihood==1;
            % compute full Gaussian  probability, necessary if the variance
            % is changing
            % TMH:UPDATE
            k = length(dd).*log(((data{id}.d_std.*sqrt(2*pi)).^(-1)));
            logL(id)=k+logL(id);
        end
        
        L(id)=exp(logL(id));
    else
        % CORRELATED GAUSSIAN
        
        N=length(data{id}.d_obs);
        
        %% CHECK OF ONE SHOULD REMOVE ICD
        if isfield(data{id},'recomputeCD');
            if (data{id}.recomputeCD==1)
                try data{id}=rmfield(data{id},'CD');end
                try data{id}=rmfield(data{id},'iCD');end
            end
        end
        
        %% CHECK WHETHER DATA SIZE HAS CHANGED
        try
            if ~(size(data{id}.iCD,1)==length(data{id}.i_use))
                % recompute CD
                try data{id}=rmfield(data{id},'CD');end
                try data{id}=rmfield(data{id},'iCD');end
                %disp('recomputing  iCD');;
            end
        end
        
        % MAKE SURE GAUSSIAN NOISE MODEL IS PROPERLY SET
        if (~isfield(data{id},'iCD'))||(data{id}.recomputeCD==1)
            if ~isfield(data{id},'CD')
                if ~isfield(data{id},'Cd');
                    if isfield(data{id},'d_std')
                        d_var=ones(N,1);
                        d_var=d_var(:).*data{id}.d_std(:).^2;
                        data{id}.Cd=diag(d_var);
                    elseif isfield(data{id},'d_var')
                        d_var=ones(N,1);
                        d_var=d_var(:).*data{id}.d_var(:);
                        data{id}.Cd=diag(d_var);
                    else
                        data{id}.Cd=0; % NO MEASUREMENT INC
                    end
                end
                if (size(data{id}.Cd,1)==1)
                    data{id}.Cd=eye(length(data{id}.d_obs)).*data{id}.Cd;
                end
                
                if ~isfield(data{id},'Ct')
                    % modelization error
                    data{id}.Ct=zeros(size(data{id}.Cd));
                end
                
                if ~isfield(data{id},'CD')
                    % modelization and measuremnet error
                    data{id}.CD=data{id}.Ct+data{id}.Cd;
                end
            end
        end
        
        if isfield(data{id},'dt');
            if length(data{id}.dt)==1
                if data{id}.use_log==1;
                    dd=(log(data{id}.d_obs(data{id}.i_use))-data{id}.dt)-log(d{id}(data{id}.i_use));
                else
                    dd=(data{id}.d_obs(data{id}.i_use)-data{id}.dt)-d{id}(data{id}.i_use);
                end
            else
                if data{id}.use_log==1;
                    %dd=(log(data{id}.d_obs(data{id}.i_use))-data{id}.dt(data{id}.i_use))-log(d{id}(data{id}.i_use));
                    dd=(log(data{id}.d_obs(data{id}.i_use))-data{id}.dt(data{id}.i_use))-log(d{id});
                else
                    %dd=(data{id}.d_obs(data{id}.i_use)-data{id}.dt(data{id}.i_use))-d{id}(data{id}.i_use);
                    dd=(data{id}.d_obs(data{id}.i_use)-data{id}.dt(data{id}.i_use))-d{id};                    
                end
            end
        else
            if data{id}.use_log==1;
                dd=log(data{id}.d_obs(data{id}.i_use))-log(d{id}(data{id}.i_use));                                
            else
                dd=data{id}.d_obs(data{id}.i_use)-d{id}(data{id}.i_use);
            end
        end
        
        if isfield(data{id},'n_score');
            % Perform forward normal score transformation of residuals to
            % normal score space           
            dd_ns = nscore_mul(dd,data{id}.n_score(data{id}.i_use));               
            dd = dd_ns;            
        end
        
        % Only compute iCD if it is computed only once (i.e.
        % data{id}.recomputeCD==0)
        if (~isfield(data{id},'iCD'))&&(data{id}.recomputeCD==0)
            %data{id}.iCD=inv(data{id}.CD);
            data{id}.iCD=inv(data{id}.CD(data{id}.i_use,data{id}.i_use));
            
            % optionally clean up data structure to save memory
            %data{1}=rmfield(data{1},'Cd');
            %data{1}=rmfield(data{1},'Ct');
            %data{1}=rmfield(data{1},'CD');
        
        end
        
        % compute logdet(CD) if it does not exist, and if recomputeCD=1;
        if (data{id}.full_likelihood)
            if (~isfield(data{id},'logdet'))||(data{id}.recomputeCD==1)
                data{id}.logdet = logdet(data{id}.CD(data{id}.i_use,data{id}.i_use));
            end
        end
        
    end %%%%%%%%%%%%%%%%%
    
    if strcmp(data{id}.noise_model,'gaussian')&&(data{id}.noise_uncorr==0)
        nknown=length(data{id}.i_use);
        
        
        if data{id}.full_likelihood==1
            f1 = -(nknown/2)*log(2*pi);
            f2 = -0.5*data{id}.logdet;
            
            if isinf(f1);
                %% this os pretty bad if CD changes !! Because then the determinant also changes..
                sippi_verbose(sprintf('%s : Full likelihood cannot be computed !',mfilename),-1)
                sippi_verbose(sprintf('%s : --> ignoring determinant !',mfilename),-1)
                f1=-f2;
            end;
            
            if data{id}.recomputeCD==1
                try
                    f3 = -.5*dd'*(data{id}.CD(data{id}.i_use,data{id}.i_use)\dd);
                catch
                    keyboard
                end
                % disp(sprintf('f3=%g, logdet=%g',f3,data{id}.logdet))
            else
                f3 =  -.5 * dd'*data{id}.iCD*dd;
            end
            %10,keyboard
            logL(id) = f1 +f2 +f3;
            
        else
            %20,keyboard
            if data{id}.recomputeCD==1
                f3 = -.5*dd'*(data{id}.CD(data{id}.i_use,data{id}.i_use)\dd);
                %f3 = -.5*dd'*(data{id}.CD\dd);
            else
                f3 =  -.5 * dd'*data{id}.iCD*dd;
            end
            %disp('===')
            %f3
            %dd(1:4)
            %data{id}.iCD(1:4,1:4)
            %data{id}.CD(1:4,1:4)
            %disp('===')
            logL(id) = f3;
            
        end
        
    elseif strcmp(data{id}.noise_model,'laplace')
        if ~isfield(data{id},'sigma');
            data{id}.sigma=sqrt(diag(data{id}.CD(data{id}.i_use,data{id}.i_use)))';
        end
        logL(id) = -.5 * sum(abs(dd(:))./data{id}.sigma(:));
    elseif (strcmp(data{id}.noise_model,'generalized_gaussian'))
        if isfield(data{id},'var');
            data{id}.sigma=sqrt(data{id}.var);
        end
        if ~isfield(data{id},'sigma');
            data{id}.sigma=sqrt(diag(data{id}.CD(data{id}.i_use,data{id}.i_use)))';
        end
        if ~isfield(data{id},'norm');
            data{id}.norm=2;
        end
        logL(id) = sum((abs(dd(:)).^data{id}.norm)./(data{id}.sigma(:).^(data{id}.norm)));
        logL(id) = logL .* (-1./data{id}.norm );
    elseif (strcmp(data{id}.noise_model,'gaussian'))
        % ALLEADY DONE
    else
        sippi_verbose(sprintf('%s : noise model ''%s'' is not supported',mfilename,data{id}.noise_model),-10);
        %keyboard
    end
    
end

if nargout>1
    L=logL;
end

logL=sum(logL);
