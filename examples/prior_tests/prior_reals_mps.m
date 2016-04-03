% prior_reals_mps: Sampling MPS type prior models
rng('default');rng(1);
clear all,close all;
ip=1; 
prior{ip}.type='mps'; % MPS type
prior{ip}.method='mps_snesim_tree'; 
prior{ip}.x=[0:.1:10]; % X array 
prior{ip}.y=[0:.1:20]; % Y array 
TI=channels;
prior{1}.ti=TI;
prior{1}.MPS.template_size=[9 9 1];

%% UNCONDITIONAL
figure(10);clf
for i=1:5;
    [m,prior]=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    axis image
    %xlabel('X')
    %ylabel('Y')
    title(sprintf('t=%3.2g',prior{1}.time))
end
colormap(sippi_colormap(1));
colorbar_shift;
print_mul(sprintf('prior_reals_%ss',prior{1}.type));
s=suptitle(sprintf('Independent realizations from a %s (%s) type prior',upper(prior{1}.type),prior{1}.method))
set(s,'interpreter','none')

%% CONDITIONAL TO HARD DATA

%% CONDITIONAL TO SOFT DATA

%% CONDITIONAL TO HARD AND SOFT DATA

