% sippi_prior_set_steplength Set step length for Metropolis sampler in SIPPI
%
% Call
%   prior=sippi_prior_set_steplength(prior,mcmc,im);
%
function  prior=sippi_prior_set_steplength(prior,mcmc,im);
if nargin<3
    im=1;
end

if ~isfield(mcmc,'i')
    mcmc.i=length(mcmc.acc(im,:));
end

i_perturb=find(mcmc.perturb(im,1:mcmc.i));
if isempty(i_perturb);
    return
end
P_current = sippi_compute_acceptance_rate(mcmc.acc(im,i_perturb),prior{im}.seq_gibbs.n_update_history);
if P_current==0, P_current=1./prior{im}.seq_gibbs.n_update_history;end
step_old=prior{im}.seq_gibbs.step;
step=sippi_adjust_step_size(step_old,P_current,prior{im}.seq_gibbs.P_target);
if (step>prior{im}.seq_gibbs.step_max), step=prior{im}.seq_gibbs.step_max;end
if (step<prior{im}.seq_gibbs.step_min), step=prior{im}.seq_gibbs.step_min;end
prior{im}.seq_gibbs.step=step;
txt_1=sprintf(' %g ',step_old);
txt_2=sprintf(' %g ',prior{im}.seq_gibbs.step);
sippi_verbose(sprintf('%s im=%02d adjusting SeqGibbs step size [%s]->[%s]. P_current=%g, P_target=%g',mfilename,im,txt_1,txt_2,P_current,prior{im}.seq_gibbs.P_target),1)
