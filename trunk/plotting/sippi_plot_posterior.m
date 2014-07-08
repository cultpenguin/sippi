function sippi_plot_posterior(fname,im_arr,prior,options,n_reals);
% sippi_plot_posterior Plot statistics from posterior sample
%
% Call :
%    sippi_plot_posterior(fname,im_arr,prior,options,n_reals);
%
% See also sippi_plot_prior
%

if nargin==0;
    [f1,fname]=fileparts(pwd);
end

options.null='';
pl_logL=0;
pl_sample=1;
pl_2d_marg=0;
pl_data=0;
pl_movie=0;

cwd=pwd;


%% DATA
if isstr(fname)
    try
        cd(fname);
        load([fname,'.mat']);
    catch
        load([fname,'.mat']);
    end

    if exist('mcmc','var')
        options.mcmc=mcmc;
    end

    
else
    data=fname;
    fname='lsq';
end


% SET DFAULT PLOTTING SETTINGS

options=sippi_plot_defaults(options);

plotdir=pwd;
try
    fname=options.txt;
end

%% PERHAPS SET THE OPTIONS IN A SEPERATE MFILE
prior=sippi_prior_init(prior);

%% logL curve
if pl_logL==1;
    sippi_plot_posterior_loglikelihood(options,prior,data,mcmc);
end

if pl_sample==1;
    sippi_plot_posterior_sample(options,prior,data,forward);
end
%% 2D POSTERIOR MARGINALS.
if (length(prior)<2); pl_2d_marg=0;end
if (pl_2d_marg==1),
    sippi_plot_posterior_2d_marg(options,prior,data);    
end

%% PLOT DATA ASSOCIATED TO REALS
if pl_data==1,
    sippi_plot_posterior_data(options,prior,data,forward);
end

%% PLOT PRIOR AND POSTERIOR MOVIE
if pl_movie==1,
    sippi_plot_movie(options.txt);
end


%%
cd(cwd);



