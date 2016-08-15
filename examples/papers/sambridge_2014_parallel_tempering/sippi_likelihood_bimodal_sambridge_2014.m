% sippi_likelihood_bimodal_sambridge_2013
%
% Likelihood accroding to Eqn. (13) in 
% Sambridge, 2014. A Parallel Tempering algorithm for probabilistic
% sampling and multimodal optimizationexample from Sambridge (2013). 
% doi: 10.1093/gji/ggt342
%
% See also: sippi_likelihood_bimodal_sambridge_2014
% 
function [logL,logL_all,data]=sippi_likelihood_bimodal_sambridge_2013(d,data,id_array)

if  nargin<3
    id_array=1:length(d);
end

for id=id_array;
    
    logL_all(id)=log(2.^(-d{id})+2.^(-(100-d{id})));
    logL_all(find(d{id}<0))=-Inf;
    logL_all(find(d{id}>100))=-Inf;
    
end

logL=sum(logL_all);