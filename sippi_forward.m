% sippi_forward Simple forward wrapper for SIPPI
%
% Assumes that the actual forward solver has been defined by
% forward.forward_function
%
% Call:
%   [d,forward,prior,data]=sippi_forward(m,forward)
%
% Optional: 
%   [d,forward,prior,data]=sippi_forward(m,forward,prior)
%   [d,forward,prior,data]=sippi_forward(m,forward,prior,data)
%   [d,forward,prior,data]=sippi_forward(m,forward,prior,data,options)
%
function [d,forward,prior,data,options]=sippi_forward(m,forward,prior,data,options)

if nargin<4;    data{1}.null='';end

% make sure to initilize the prior if it has not allready been done
% TMH: can this be ignored?
if exist('prior','var')
for ip=1:length(prior)
    if ~isfield(prior{ip},'init');
        prior=sippi_prior_init(prior);
    end
end
end

if isfield(forward,'forward_function');
    
    if nargin==2;
        [d,forward]=feval(forward.forward_function,m,forward);
    elseif nargin==3
        [d,forward,prior]=feval(forward.forward_function,m,forward,prior);
    elseif nargin==4

        [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior,data);
    elseif nargin==5
        [d,forward,prior,data,options]=feval(forward.forward_function,m,forward,prior,data,options);
    end
else
    disp(sprintf('%s : No forward_function specified in ''forward'' structure',mfilename))
    d=[];
    
end

%% Optionally adjust the noise model

if nargin>=4
    for im=1:length(prior);
        
        if  ~isfield(prior{im},'is_perturbed');
            prior{im}.is_perturbed=0;
        end
        
        % adjust uncorrelated noise
        if ( (strcmp(lower(prior{im}.name),'d_std')||strcmp(lower(prior{im}.name),'d_var')) && (prior{im}.is_perturbed==1));
            
            if  ~isfield(prior{im},'perturb_noise');
                prior{im}.perturb_noise=0;
            end
            if  ~isfield(prior{im},'perturb_noise_i_data');
                prior{im}.perturb_noise_i_data=1;
            end
            id=prior{im}.perturb_noise_i_data;
            
            % UPDATE UNCORREALTED NOISE MODEL
            if prior{im}.perturb_noise==1;
                if isfield(data{id},'Cd');
                    data{id}=rmfield(data{id},'Cd');
                end
                
                if strcmp(lower(prior{im}.name),'d_std') 
                    data{id}.d_std=m{im};
                    sippi_verbose(sprintf('%s: updating noise on data to d_std=%6.2g',mfilename,data{id}.d_std),2);                
                else
                    data{id}.d_var=m{im};
                    sippi_verbose(sprintf('%s: updating noise on data to d_var=%6.2g',mfilename,data{id}.d_var),2);                
                end
                data{id}.recomputeCD=1;
                data{id}.full_likelihood=1;
                
            end
        end
    
        
        % adjust uncorrelated noise
        
        if  (strcmp(lower(prior{im}.name),'ct') && (prior{im}.is_perturbed==1));
            if  ~isfield(prior{im},'perturb_noise');
                prior{im}.perturb_noise=0;
            end
            if  ~isfield(prior{im},'perturb_noise_i_data');
                prior{im}.perturb_noise_i_data=1;
            end
            id=prior{im}.perturb_noise_i_data;
            
            if prior{im}.perturb_noise==1;
                if isfield(data{id},'Ct');
                    data{id}=rmfield(data{id},'Ct');
                end
                if  ~isfield(prior{im},'Ct_base');
                    prior{im}.Ct_base=eye(length(data{1}.d_obs));
                end
                
                d_fac=m{im};
                data{id}.Ct=d_fac.*prior{im}.Ct_base;
                
                data{id}.recomputeCD=1;
                data{id}.full_likelihood=1;
                
            end
            
        end
        
        
    end
end    
    

