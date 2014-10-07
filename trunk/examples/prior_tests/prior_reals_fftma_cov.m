% prior_reals_fftma_cov Sampling an FFT-MA type prior model with varying covariance properties


clear all;
try;close(12);end
try;close(13);end

im=0;

im=im+1; 
prior{im}.type='FFTMA';
prior{im}.x=[0:.1:10]; % X array 
prior{im}.y=[0:.1:20]; % Y array 
prior{im}.m0=10;
prior{im}.Cm='1 Sph(10,90,.25)';

im=im+1;
prior{im}.type='uniform';
prior{im}.name='range_1';
prior{im}.min=2;
prior{im}.max=14;
prior{im}.prior_master=1;

im=im+1;
prior{im}.type='gaussian';
prior{im}.name='ang_1';
prior{im}.m0=90;
prior{im}.std=10;
prior{im}.prior_master=1;


%%
randn('seed',4);
figure(12);clf
for i=1:5;
    m=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    %caxis([8 12])
    axis image
    %xlabel('X')
    %ylabel('Y')
end
colorbar_shift;
colormap(sippi_colormap(1));
print_mul('prior_reals_fftma_cov')
suptitle(sprintf('FFT-MA with varying covariance properties'))


%% SEQ GIBBS
% 
randn('seed',4);rand('seed',4);
figure(13);clf
prior{1}.seq_gibbs.step=.5;
prior{2}.seq_gibbs.step=.25;
prior{3}.seq_gibbs.step=0;;
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
colormap(sippi_colormap(1));
print_mul('prior_reals_fftma_cov_seqgibbs')
suptitle(sprintf('FFT-MA with varying covariance properties\n Sequential Gibbs sampling'))


%% SEQ GIBBS MOVIE
% 
fclose all;
close all
randn('seed',4);rand('seed',4);
fig=figure(13);clf
prior{1}.seq_gibbs.step=.15;
prior{2}.seq_gibbs.step=.15;
prior{3}.seq_gibbs.step=0.05;;
prior{3}.perturb=1;
[m,prior]=sippi_prior(prior);
aviobj = avifile('prior_reals_fftma_cov_seqgibbs.avi')
for i=1:250;
    [m,prior]=sippi_prior(prior,m);
    subplot(1,1,1);
    imagesc(prior{3}.x,prior{3}.y,m{3});
    caxis([8 12])
    axis image
    %colorbar_shift;
    colormap(sippi_colormap(1));
    drawnow;
    F = getframe(fig);
    aviobj = addframe(aviobj,F);
     %xlabel('X')
    %ylabel('Y')
end
aviobj = close(aviobj);
 
%print_mul('prior_reals_fftma_cov_seqgibbs')
%suptitle(sprintf('FFT-MA with varying covariance properties\n Sequential Gibbs sampling'))


