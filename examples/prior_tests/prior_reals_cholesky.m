% prior_reals_cholesky Sampling a CHOLESKT type prior model
clear all,close all;

%% Define a 2D Cholesky type a priori model
im=1;
prior{im}.type='cholesky';
dx=0.2;
prior{im}.x=[0:dx:10]; % X array
prior{im}.y=[0:dx:20]; % Y array
prior{im}.m0=10;
prior{im}.Cm='1 Sph(10,90,.25)';

randn('seed',1);
%%  Plot sample of the prior model
figure(10);clf
for i=1:5;
    [m,prior]=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis([8 12])
    axis image
end
colormap(sippi_colormap(1));
colorbar_shift;
print_mul('prior_reals_cholesky')
suptitle(sprintf('CHOLESKY, Cm=%s',prior{im}.Cm))

%% Setup target distribution
N=10000;
prob_chan=0.5;
d1=randn(1,ceil(N*(1-prob_chan)))*.5+8.5;
d2=randn(1,ceil(N*(prob_chan)))*.5+11.5;
d_target=[d1(:);d2(:)];
% UPDATE PRIOR STRUCTURE TO USE TARGET DISTRIBUTION
prior{im}.d_target=d_target;
% using d_target, with no trend:
prior{im}.m0=0;

[m,prior]=sippi_prior(prior);

%%  Plot sample of the prior model, using target distribution
randn('seed',1);
figure(11);clf
for i=1:5;
    m=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    caxis([8 12])
    axis image
    %xlabel('X')
    %ylabel('Y')
end
colormap(sippi_colormap(1));
suptitle(sprintf('CHOLESKY using target histogram, Cm=%s',prior{im}.Cm));
colorbar_shift;
print_mul('prior_reals_cholesky_nscore')

%% Show movie of a random walk in the prior (using sequential Gibbs sampling)
save_movie=0;
for j=1:2
    randn('seed',1);
    figure(11+j);clf
    if save_movie==1;set(gcf,'units','normalized','outerposition',[0 0 1 1]);end
    prior{1}.seq_gibbs.step=0.02;
    
    if (j==1);
        if (save_movie==1)
            fname=[mfilename,'_target','.mp4'];
            vidObj = VideoWriter(fname,'MPEG-4');
        end
    else
        prior{1}=rmfield(prior{1},'d_target');
        prior{1}=rmfield(prior{1},'o_nscore');
        prior{1}.m0=10;
        if save_movie==1
            fname=[mfilename,'.mp4'];
            vidObj = VideoWriter(fname,'MPEG-4');
        end
    end
    if save_movie==1;open(vidObj);end
    for i=1:200;
        [m,prior]=sippi_prior(prior,m);
        imagesc(prior{1}.x,prior{1}.y,m{1});
        colormap(sippi_colormap(1));
        caxis([8 12])
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
end
