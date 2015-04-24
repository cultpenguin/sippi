% sippi_anneal_adjust_noise : Adjust noise level in annealing schedul
%
% Call: 
%    [data_adjust,mcmc]=sippi_anneal_adjust_noise(data,i,mcmc,prior);
%
% See also: sippi_metropolis, sippi_anneal_factor
%
function [data,mcmc]=sippi_anneal_adjust_noise(data_org,i,mcmc,prior);

data=data_org;
%% GET NOISE SCALING FACTOR
if nargin==4;
    [fac,mcmc]=sippi_anneal_factor(mcmc,i,prior);
elseif nargin==2
    fac=i;
else
    [fac,mcmc]=sippi_anneal_factor(mcmc,i);
end

%% APPLY NOISE AMPLIFICATION
for id=1:length(data_org)
    if isfield(data{id},'d_std');
        data{id}.d_std=fac.*data_org{id}.d_std;
        %disp(sprintf('%s : i=%05d D(%d), d_std_new=%5g (%5g)',mfilename,i,id,data{id}.d_std,data_org{id}.d_std))
    end
    
    if isfield(data{id},'Cd');
        data{id}.Cd=(fac.*sqrt(data_org{id}.Cd)).^2;
%        d_std=fac.*sqrt(diag(data_org{id}.Cd));;
%        for j=1:size(data{id}.Cd,1);
%            data{id}.Cd(j,j)=d_std(j).^2;
%        end
    end
    
    if isfield(data{id},'CD');
        data{id}.CD=(fac.*sqrt(data_org{id}.CD)).^2;
    end
    if isfield(data{id},'iCD');
        data{id}=rmfield(data{id},'iCD');
    end
    
    %% NOTE: WE PROBABLY NEED TO MAKE SURE THIS NOISE ADJUSTMENT FOR FOR BOTH
    % CORRELATED AND UNCORRELATED ERRORS, d_std, Cd, Ct, CD !!
end



