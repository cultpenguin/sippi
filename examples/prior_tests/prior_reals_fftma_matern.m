% prior_reals_fftma_matern Sampling an FFT-MA type prior model with varying
% covariance properties using a Matern type covariance model


clear all;
try;close(12);end
try;close(13);end

%% Setup FFT-MA type a priori model with uncertainty covariance model parameters and Matern type covariance model
im=0;

im=im+1; 
prior{im}.type='FFTMA';
prior{im}.x=[0:.1:10]; % X array 
prior{im}.y=[0:.1:20]; % Y array 
prior{im}.m0=10;
prior{im}.Va='1 Mat(5,90,.5)';
prior{im}.fftma_options.constant_C=0;


im=im+1;
prior{im}.type='gaussian';
prior{im}.name='range_1';
prior{im}.m0=1.5;
prior{im}.std=.5;
prior{im}.norm=80;
prior{im}.prior_master=1;

im=im+1;
prior{im}.type='gaussian';
prior{im}.name='ang_1';
prior{im}.m0=90;
prior{im}.std=40;
prior{im}.norm=80;
prior{im}.prior_master=1;


im=im+1;
prior{im}.type='gaussian';
prior{im}.name='nu'; % parameter controling Matern type covariance model
prior{im}.m0=1;
prior{im}.std=.9;
prior{im}.norm=80;
prior{im}.prior_master=1;

prior=sippi_prior_init(prior);

%% Plot sample of the prior model
randn('seed',4);
figure(12);clf
for i=1:5;
    m=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis([8 12])
    axis image
end
colorbar_shift;
colormap(sippi_colormap(1));
print_mul('prior_reals_fftma_matern')
suptitle(sprintf('FFT-MA with varying Matern covariance properties'))


%% Plot 5 models generated using seqential Gibbs sampling
randn('seed',4);rand('seed',4);
figure(13);clf
prior{1}.seq_gibbs.step=0;
prior{2}.seq_gibbs.step=.1;
prior{3}.seq_gibbs.step=.1;;
prior{4}.seq_gibbs.step=.4;;
prior{1}.perturb=1;
[m,prior]=sippi_prior(prior);
for i=1:5;
    [m,prior]=sippi_prior(prior,m);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis([8 12])
    axis image
    %xlabel('X')
    %ylabel('Y')
end
colorbar_shift;
colormap(sippi_colormap(1));
print_mul('prior_reals_fftma_matern_seqgibbs')
suptitle(sprintf('FFT-MA with varying Matern type covariance properties\n Sequential Gibbs sampling'))


%% Show movie of a random walk in the prior (using sequential Gibbs sampling)
save_movie=0;

fclose all;
close all
randn('seed',4);rand('seed',4);
fig=figure(13);clf
prior{2}.seq_gibbs.step=.05;
prior{3}.seq_gibbs.step=.05;
prior{4}.seq_gibbs.step=.1;.025;
prior{1}.seq_gibbs.step=.01;0.01;;
prior{1}.perturb=1;
[m,prior]=sippi_prior(prior);
if (save_movie==1);
  fname=[mfilename,'.mp4'];
  vidObj = VideoWriter(fname,'MPEG-4');
  open(vidObj);
end

for i=1:250;
    [m,prior]=sippi_prior(prior,m);
    subplot(1,1,1);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis([8 12])
    axis image
    %colorbar_shift;
    colormap(sippi_colormap(1));
    drawnow;
    if save_movie==1
        % Write each frame to the file.
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
end
if save_movie==1;close(vidObj);end
