% sippi_likelihood_multinomial: Compute likelihood using the multinomial function
%
% Call
%   [logL,L,data]=sippi_likelihood_multinomial(d,data);
%
% See also sippi_likelihood
%
function [logL,logL_all,data]=sippi_likelihood_multinom(d,data,id_array);
if  nargin<3
    id_array=1:length(d);
end

for id=id_array
    
    logL_all(id)=0;

    if ~isfield(data{id},'nprior');
        data{id}.nprior=0;
    end
    
    logL_all(id)=multinomial(d{id},data{id}.d_obs,data{id}.nprior);
end
logL=sum(logL_all);

