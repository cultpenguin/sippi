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
    
    logL_all(id)=multinomial(d{id},data{id}.d_obs);
end
logL=sum(logL_all);

