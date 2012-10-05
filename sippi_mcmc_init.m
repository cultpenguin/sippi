%
function options=sippi_mcmc_init(options,prior);

options.mcmc.null='';
if ~isfield(options.mcmc,'nite');options.mcmc.nite=30000;end
if ~isfield(options.mcmc,'i_sample');options.mcmc.i_sample=500;end
if ~isfield(options.mcmc,'i_plot');options.mcmc.i_plot=50;end

if ~isfield(options.mcmc,'accept_only_improvements');options.mcmc.accept_only_improvements=0;end
if ~isfield(options.mcmc,'accept_all');options.mcmc.accept_all=0;end

if ~isfield(options.mcmc,'store_all');options.mcmc.store_all=0;end


% pertubation strategy
if nargin>1
    if ~isfield(options.mcmc,'pert_strategy');
        options.mcmc.pert_strategy.i_pert=1:1:length(prior);
    end
    if ~isfield(options.mcmc.pert_strategy,'i_pert');
        options.mcmc.pert_strategy.i_pert=1:1:length(prior);
    end
    if ~isfield(options.mcmc.pert_strategy,'i_pert_freq');
        options.mcmc.pert_strategy.i_pert_freq=ones(1,length(prior))./length(prior);
    end
    
    if ~isfield(options.mcmc.pert_strategy,'perturb_all');
        options.mcmc.pert_strategy.perturb_all=0;
    end
    
end
