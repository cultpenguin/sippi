% prior_reals_fftma_cov Sampling an FFT-MA type prior model with varying covariance properties


clear all
im=1;
prior{im}.type='gaussian';
prior{im}.name='range_1';
prior{im}.m0=10;
prior{im}.std=4;
prior{im}.norm=80;
prior{im}.prior_master=3;

im=2;
prior{im}.type='gaussian';
prior{im}.name='ang_1';
prior{im}.m0=90;
prior{im}.std=50;
prior{im}.norm=80;
prior{im}.prior_master=3;

im=3; 
prior{im}.type='FFTMA';
prior{im}.x=[0:.1:10]; % X array 
prior{im}.y=[0:.1:20]; % Y array 
prior{im}.m0=10;
prior{im}.Va='1 Sph(10,90,.25)';
prior{im}.fftma_options.constant_C=0;

prior=sippi_prior_init(prior);

%%
randn('seed',4);
figure(12);clf
for i=1:5;
    m=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{3}.x,prior{3}.y,m{3});
    caxis([8 12])
    axis image
    %xlabel('X')
    %ylabel('Y')
end
colorbar_shift;
try;load cmap;colormap(cmap);end
print_mul('prior_reals_fftma_cov')
suptitle(sprintf('FFT-MA with varying covariance properties'))


%% SEQ GIBBS
% 
randn('seed',4);rand('seed',4);
figure(13);clf
prior{1}.seq_gibbs.step=.5;
prior{2}.seq_gibbs.step=.25;
prior{3}.seq_gibbs.step=0;
prior{3}.perturb=1;
[m,prior]=sippi_prior(prior);
for i=1:5;
    [m,prior]=sippi_prior(prior,m);
    subplot(1,5,i);
    imagesc(prior{3}.x,prior{3}.y,m{3});
    caxis([8 12])
    axis image
    %xlabel('X')
    %ylabel('Y')
end
colorbar_shift;
try;load cmap;colormap(cmap);end
print_mul('prior_reals_fftma_cov_seqgibbs')
suptitle(sprintf('FFT-MA with varying covariance properties\n Sequential Gibbs sampling'))


