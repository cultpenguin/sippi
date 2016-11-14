% prior_reals_plurigaussian Sampling an a priori model based on voronoi
% cells
clear all,close all;
rng('default')

%% Define a 2D a priori model based on 10 VORONOIS CELLS
cells_N_max=50;
ip=1;
prior{ip}.type='voronoi';
prior{ip}.x=1:.04:20;prior{1}.y=1:.04:20;
prior{ip}.cells_N=cells_N_max; % SET NUMBER OF CELLS
prior{ip}.cax=[1 prior{1}.cells_N];
[m,prior]=sippi_prior(prior);
sippi_plot_prior(prior,m);


%%  Plot sample of the prior model
figure(10);clf
for i=1:5;
    [m,prior]=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis(prior{1}.cax)
    axis image
end
colormap(sippi_colormap(1));
colorbar_shift;
drawnow;
print_mul(mfilename)
%suptitle(sprintf('voronois'))


%% Define a 2D a priori model based on 10-30 VORONOI CELLS

ip=ip+1;
prior{ip}.type='uniform';
prior{ip}.name='cells_N';
prior{ip}.min=5;
prior{ip}.max=50;
prior{ip}.prior_master=1;
prior{1}.cax=[1 50];

ip=ip+1;
prior{ip}.type='uniform';
prior{ip}.name='cells_x';
prior{ip}.x=[1:cells_N_max];
prior{ip}.min=min(prior{1}.x);
prior{ip}.max=max(prior{1}.x);
prior{ip}.prior_master=1;

ip=ip+1;
prior{ip}.type='uniform';
prior{ip}.name='cells_y';
prior{ip}.x=[1:cells_N_max];
prior{ip}.min=min(prior{1}.y);
prior{ip}.max=max(prior{1}.y);
prior{ip}.prior_master=1;

ip=ip+1;
prior{ip}.type='uniform';
prior{ip}.name='cells_value';
prior{ip}.x=[1:cells_N_max];
prior{ip}.min=5;
prior{ip}.max=15;
prior{ip}.prior_master=1;
prior{1}.cax=[5 15];

rng(1);
figure(11);clf
for i=1:5;
    [m,prior]=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis(prior{1}.cax)
    axis image
    title(sprintf('cells_N=%d',prior{1}.cells_N))
end
colormap(sippi_colormap(1));
colorbar_shift;
drawnow;
print_mul([mfilename,'_cells_N'])

%% Show movie of a random walk in the prior (using sequential Gibbs sampling)
save_movie=1;
randn('seed',1);
figure(11+j);clf
if save_movie==1;set(gcf,'units','normalized','outerposition',[0 0 1 1]);end
for ip=1:length(prior)
    prior{ip}.seq_gibbs.step=0.02;
end
if (save_movie==1)
    if isunix==1
        fname=[mfilename,'_target','.avi'];
        vidObj = VideoWriter(fname,'Uncompressed AVI');
    else
        fname=[mfilename,'_target','.mp4'];
        vidObj = VideoWriter(fname,'MPEG-4');
    end
end

if save_movie==1;open(vidObj);end
for i=1:500;
    [m,prior]=sippi_prior(prior,m);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    colormap(sippi_colormap(1));
    caxis(prior{1}.cax)
    axis image
    axis tight
    set(gca,'nextplot','replacechildren');
    drawnow;
    if save_movie==1
        % Write each frame to the file.
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
end
if save_movie==1;close(vidObj);end
