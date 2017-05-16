% prior_reals_mps: Sampling MIXSIM type prior models
rng('default');rng(1);
clear all,close all;


ip=1; 
prior{ip}.type='mixsim'; % MPS type
%prior{ip}.x=[0:.5:10]; % X array 
%prior{ip}.y=[0:.5:20]; % Y array 
prior{ip}.x=[1:1:20]; % X array 
prior{ip}.y=[1:1:40]; % Y array 

prior{ip}.x=[1:1:60]; % X array  % X array (CURRENTLY ONLY FOR INTEGEGER X STARTING FROM 1, WITH DX=1)
prior{ip}.y=[1:1:20]; % Y array  % Y array (CURRENTLY ONLY FOR INTEGEGER Y STARTING FROM 1, WITH DY=1)



%% set TI
% Strebelles channels
TI=channels;
TI=TI(3:3:100,1:1:10)+1;
% 1D logs
%clear TI;load Logs1D.mat;for mm=1:size(Logs1D,2);TI{mm}=Logs1D(1:end,mm);end
    
prior{ip}.TI=TI;





%%
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
    drawnow;
end
caxis([1 prior{1}.options.sV])
colormap(sippi_colormap(1));
colorbar_shift;
print_mul(sprintf('prior_reals_%ss',prior{1}.type));
s=suptitle(sprintf('Independent realizations from a %s type prior',prior{1}.type))
set(s,'interpreter','none')


%% SEQUENTIAL GIBBS

[m2,prior]=sippi_prior(prior,m);

%% SEQUENTIAL GIBBS
figure(11);clf
prior{1}.seq_gibbs.type=[2];
prior{1}.seq_gibbs.step=[.9];

%prior{1}.seq_gibbs.type=[1];
%prior{1}.seq_gibbs.step=[10];

for i=1:5;
    [m,prior]=sippi_prior(prior,m);
    figure(11);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    axis image
    colormap(sippi_colormap(1));
    caxis([1 prior{1}.options.sV])
    axis tight
    set(gca,'nextplot','replacechildren');
    drawnow;
end
print_mul(sprintf('prior_reals_%s_%d',prior{1}.type,prior{1}.seq_gibbs.type));
s=suptitle(sprintf('Independent realizations from a %s type prior (%d)',prior{1}.type,prior{1}.seq_gibbs.type))
set(s,'interpreter','none')
    
return
%% PRIOR MOVIE
fclose all;
fname=[mfilename,'_target','.mp4'];
vidObj = VideoWriter(fname,'MPEG-4');
open(vidObj);
N=200;
for i=1:N;  
    [m,prior]=sippi_prior(prior,m);
    sippi_plot_prior(prior,m);
    colormap(sippi_colormap(1));
    caxis([1 prior{1}.options.sV]);
    title(sprintf('Realizations %03d/%03d',i,N));
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
end
close(vidObj);