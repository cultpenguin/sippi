% prior_reals_snesim Sampling a SNESIM type prior model
clear all,
ip=1; 
prior{ip}.type='SNESIM'; % FORTRAN VERSION OF SNESIM
%prior{ip}.type='SNESIM_STD'; % SNESIM FROM SGeMS
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
print_mul(sprintf('prior_reals_%s',prior{1}.type));
suptitle(sprintf('Independent samples from a SNESIM type prior'))


%% 
figure(11);clf
[m,prior]=sippi_prior(prior);
for i=1:5;
    try;prior{ip}.S.XML.parameters.Nb_Multigrids_ADVANCED.value=i;end;
    try;prior{ip}.S.nmulgrids=i;end
    m=sippi_prior(prior);
	subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});axis image	
    title(sprintf('NM GRIDS = %d',i))
end
colormap(sippi_colormap(1));
colorbar_shift;
print_mul(sprintf('prior_reals_%s_nmgrid',prior{1}.type));


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
print_mul(sprintf('prior_reals_%s_rotation_scale',prior{1}.type));
suptitle(sprintf('Independant samples from a SNESIM type prior'))



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
print_mul(sprintf('prior_reals_%s_seqgibbs_type1',prior{ip}.type));
suptitle(sprintf('Sampling from a SNESIM type prior\nSequential Gibbs sampling type 1'))


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
print_mul(sprintf('prior_reals_%s_seqgibbs_type2',prior{ip}.type));
suptitle(sprintf('Sampling from a SNESIM type prior\nSequential Gibbs sampling type 2'))

sgems_clean;

