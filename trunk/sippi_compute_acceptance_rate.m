% sippi_compute_acceptance_rate Computes acceptance rate for the Metropolis sampler in SIPPI
%
% Call:
%   P_acc=sippi_compute_acceptance_rate(acc,n_update_history);
%
function P_acc=sippi_compute_acceptance_rate(acc,n_update_history);
if nargin<3
    n_update_history=50;
end

i_max=length(acc);
i1=max([i_max-n_update_history,1]);
ii=i1:1:i_max;
P_acc=sum(acc(ii))./length(ii);
