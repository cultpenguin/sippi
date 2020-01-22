% sippi_mcmc_init Initialize McMC options for Metropolis and rejection sampling in SIPPI
%
% Call:
%    options=sippi_mcmc_init(options,prior);
%
function options=sippi_mcmc_init(options,prior);

options.mcmc.null='';
if ~isfield(options.mcmc,'nite');options.mcmc.nite=30000;end
if isfield(options.mcmc,'n_sample');options.mcmc.i_sample=max([1,ceil(options.mcmc.nite/options.mcmc.n_sample)]);end
if ~isfield(options.mcmc,'i_sample');options.mcmc.i_sample=100;end
if ~isfield(options.mcmc,'i_save_workspace');options.mcmc.i_save_workspace=options.mcmc.nite;end
if ~isfield(options.mcmc,'i_plot');options.mcmc.i_plot=50;end

if ~isfield(options.mcmc,'accept_only_improvements');options.mcmc.accept_only_improvements=0;end
if ~isfield(options.mcmc,'accept_all');options.mcmc.accept_all=0;end

if ~isfield(options.mcmc,'store_all');options.mcmc.store_all=0;end

% chains
if ~isfield(options,'mcmc'); options.mcmc.null='';end
if ~isfield(options.mcmc,'n_chains'); options.mcmc.n_chains=1;end

%% pertubation strategy
if nargin>1
    if ~isfield(options.mcmc,'pert_strategy');
        options.mcmc.pert_strategy.i_pert=1:1:length(prior);
    end
    if ~isfield(options.mcmc.pert_strategy,'i_pert');
        options.mcmc.pert_strategy.i_pert=1:1:length(prior);
    end
    if ~isfield(options.mcmc.pert_strategy,'i_pert_freq');
        np=length(options.mcmc.pert_strategy.i_pert);
        options.mcmc.pert_strategy.i_pert_freq=ones(1,np)./np;
    end
    
    if ~isfield(options.mcmc.pert_strategy,'perturb_all');
        options.mcmc.pert_strategy.perturb_all=0;
    end
    
end

%% Gibbs type optimization
options.mcmc.gibbs.null='';
if ~isfield(options.mcmc.gibbs,'usedim'); options.mcmc.gibbs.usedim=1;end
if ~isfield(options.mcmc.gibbs,'i_gibbs'); options.mcmc.gibbs.i_gibbs=1e+9;end
if ~isfield(options.mcmc.gibbs,'Nm');
    if (options.mcmc.gibbs.usedim==1);
        options.mcmc.gibbs.N_bins=41;
    else
        options.mcmc.gibbs.N_bins=300;
    end
end

%% set Temperature
if ~isfield(options.mcmc,'T');
    if ~isfield(options.mcmc,'T_min');options.mcmc.T_min=1;end
    if ~isfield(options.mcmc,'T_max');options.mcmc.T_max=options.mcmc.n_chains;end
    
    if options.mcmc.n_chains==1;
        options.mcmc.T=options.mcmc.T_min;
    else
        options.mcmc.T=linspace(options.mcmc.T_min,options.mcmc.T_max,NC);
    end
end

%% parallel tempering frequency jump
if ~isfield(options.mcmc,'chain_frequency_jump');
    options.mcmc.chain_frequency_jump=0.1;
end

%% Parallel, ordering of chains
if ~isfield(options.mcmc,'order_chains_frequency');
    options.mcmc.order_chains_frequency=0;
end
%% CHECK FOR ANNEALING
if isfield(options.mcmc,'anneal');
    options.mcmc.do_anneal=1;
    options.mcmc.T_fac=1;
else
    options.mcmc.do_anneal=0;
    options.mcmc.T_fac=1;
end


