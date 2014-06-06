% prior_reals_snesim Sampling a SNESIM type prior model
clear all,
ip=1; 
prior{ip}.type='SNESIM';
prior{ip}.x=[0:.1:10]; % X array 
prior{ip}.y=[0:.1:20]; % Y array 

%% get Training image from file
% f_ti='3194256336_e75579cd07_m.sgems';
% %f_ti='2485773707_997e223aab_n.sgems';
% O=sgems_read(f_ti);
% ti=O.D(:,:,:,3);figure(1);;imagesc(ti);axis image;colorbar;title('ti')
% prior{ip}.ti=ti;
% %
% prior{ip}.index_values=[0 1 2];
% prior{ip}.m_values=[8 12 10];


randn('seed',1);rand('seed',1);

%% UNCONDITIONAL NO SCALING
figure(10);clf
for i=1:5;
    m=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    axis image
    %xlabel('X')
    %ylabel('Y')
end
colormap(sippi_colormap(1));
colorbar_shift;
print_mul('prior_reals_snesim');
suptitle(sprintf('Independant samples from a SNESIM type prior'))


%% 
figure(11);clf
[m,prior]=sippi_prior(prior);
for i=1:5;
    prior{ip}.S.XML.parameters.Nb_Multigrids_ADVANCED.value=i;
    m=sippi_prior(prior);
	subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});axis image	
    title(sprintf('NM GRIDS = %d',i))
end
colormap(sippi_colormap(1));
colorbar_shift;
print_mul('prior_reals_snesim_nmgrid');


%% UNCONDITIONAL WITH SCALING
figure(12);clf
prior{ip}.scaling=.7;
prior{ip}.rotation=30;
for i=1:5;
    m=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    axis image
    %xlabel('X')
    %ylabel('Y')
end
colormap(sippi_colormap(1));
colorbar_shift;
print_mul('prior_reals_snesim_rotation_scale');
suptitle(sprintf('Independant samples from a SNESIM type prior'))



return
%% SEQ SIM 1
randn('seed',1);rand('seed',1);
prior{ip}.seq_gibbs.type=1;
prior{ip}.seq_gibbs.step=[20];
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
colormap(sippi_colormap(1));
colorbar_shift;
print_mul('prior_reals_snesim_seqgibbs_type1');
suptitle(sprintf('Sampling from a SNESIM type prior\nSequential Gibbs sampling'))


%% SEQ SIM 2
randn('seed',1);rand('seed',1);
prior{ip}.seq_gibbs.type=2;
prior{ip}.seq_gibbs.step=.97;
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
colormap(sippi_colormap(1));

colorbar_shift;
print_mul('prior_reals_snesim_seqgibbs_type2');
suptitle(sprintf('Sampling from a SNESIM type prior\nSequential Gibbs sampling'))

sgems_clean;

