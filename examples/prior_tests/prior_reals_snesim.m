% prior_reals_snesim Sampling a SNESIM type prior model
clear all,
im=1; 
prior{im}.type='SNESIM';
prior{im}.x=[0:.1:10]; % X array 
prior{im}.y=[0:.1:20]; % Y array 

%% get Training image from file
% f_ti='3194256336_e75579cd07_m.sgems';
% %f_ti='2485773707_997e223aab_n.sgems';
% O=sgems_read(f_ti);
% ti=O.D(:,:,:,3);figure(1);;imagesc(ti);axis image;colorbar;title('ti')
% prior{im}.ti=ti;
% %
% prior{im}.index_values=[0 1 2];
% prior{im}.m_values=[8 12 10];


prior{im}.S=sgems_get_par('snesim_std');
prior=sippi_prior_init(prior);

%prior{im}.scaling=.7;
prior{im}.rotation=30;


randn('seed',1);rand('seed',1);

O=sgems_grid(prior{1}.S)

[m,prior]=sippi_prior(prior);
sippi_plot_model(prior,m)
return

%% UNCONDITIONAL
figure(10);clf
for i=1:5;
    m=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    axis image
    %xlabel('X')
    %ylabel('Y')
end
try;load cmap;colormap(cmap);end
colorbar_shift;
print_mul('prior_reals_snesim');
suptitle(sprintf('Independant samples from a SNESIM type prior'))



%% SEQ SIM 1
randn('seed',1);rand('seed',1);
prior{im}.seq_gibbs.type=1;
prior{im}.seq_gibbs.step=[20];
[m,prior]=sippi_prior(prior);
[m,prior]=sippi_prior(prior,m);  
figure(11);clf
for i=1:5;
    m_old=m;
    [m,prior]=sippi_prior(prior,m);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    axis image
    %xlabel('X')
    %ylabel('Y')
end
try;load cmap;colormap(cmap);end
colorbar_shift;
print_mul('prior_reals_snesim_seqgibbs_type1');
suptitle(sprintf('Sampling from a SNESIM type prior\nSequential Gibbs sampling'))


%% SEQ SIM 2
randn('seed',1);rand('seed',1);
prior{im}.seq_gibbs.type=2;
prior{im}.seq_gibbs.step=.97;
[m,prior]=sippi_prior(prior);
[m,prior]=sippi_prior(prior,m);  
figure(12);clf
for i=1:5;
    m_old=m;
    [m,prior]=sippi_prior(prior,m);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    axis image
    %xlabel('X')
    %ylabel('Y')
end
try;load cmap;colormap(cmap);end

colorbar_shift;
print_mul('prior_reals_snesim_seqgibbs_type2');
suptitle(sprintf('Sampling from a SNESIM type prior\nSequential Gibbs sampling'))

sgems_clean

