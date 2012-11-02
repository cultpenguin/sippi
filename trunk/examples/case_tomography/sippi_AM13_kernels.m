% sippi_AM13_kernels : plot som sensitivity kernels associated to forward type

clear all;close all
D=load('AM13_data.mat');

%% THE PRIOR
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.m0=0.145;
prior{im}.Va='.0003 Sph(6)';

dx=0.125;
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];
prior{im}.cax=[0.1 0.18];
prior=sippi_prior_init(prior);
%% THE DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;


%% FORWARD MODEL
D=load('AM13_data.mat');
forward.sources=D.S;
forward.receivers=D.R;
forward.forward_function='sippi_forward_traveltime';

%% KERNELS
data{1}.i_use=[1,100,200];
m=sippi_prior(prior);

j=0;
forward.freq=.1;

j=j+1;forward.type='ray';forward.linear=1;
[d{j},forward_mul{j}]=sippi_forward(m,forward,prior,data,id,im);
K{j} = reshape(sum(forward_mul{j}.G(:,:)),prior{1}.dim(2),prior{1}.dim(1));
txt{j}='Straight ray';

j=j+1;forward.type='ray';forward.linear=0;
[d{j},forward_mul{j}]=sippi_forward(m,forward,prior,data,id,im);
K{j} = reshape(sum(forward_mul{j}.G(:,:)),prior{1}.dim(2),prior{1}.dim(1));
txt{j}='Bended ray';

j=j+1;forward.type='fat';forward.linear=1;
[d{j},forward_mul{j}]=sippi_forward(m,forward,prior,data,id,im);
K{j} = reshape(sum(forward_mul{j}.G(:,:)),prior{1}.dim(2),prior{1}.dim(1));
txt{j}='Straight fat beam';

j=j+1;forward.type='fat';forward.linear=0;
[d{j},forward_mul{j}]=sippi_forward(m,forward,prior,data,id,im);
K{j} = reshape(sum(forward_mul{j}.G(:,:)),prior{1}.dim(2),prior{1}.dim(1));
txt{j}='Bended fat beam';

j=j+1;forward.type='born';forward.linear=1;
[d{j},forward_mul{j}]=sippi_forward(m,forward,prior,data,id,im);
K{j} = reshape(sum(forward_mul{j}.G(:,:)),prior{1}.dim(2),prior{1}.dim(1));
txt{j}='Linear Born kernel';

j=j+1;forward.type='born';forward.linear=0;
[d{j},forward_mul{j}]=sippi_forward_traveltime(m,forward,prior,data,id,im);
K{j} = reshape(sum(forward_mul{j}.G(:,:)),prior{1}.dim(2),prior{1}.dim(1));
txt{j}='Non-linear Born kernel';




%% PLOT KERNELS
figure(1);set_paper;clf;
subplot(1,7,1);
imagesc(prior{1}.x,prior{1}.y,m{1});axis image
title('Reference model')

for i=1:length(K);
    figure_focus(i+1);
    imagesc(prior{1}.x,prior{1}.y,K{i});axis image
    title(txt{i})
    caxis([-1 1].*.02)
    colorbar
    print_mul(sprintf('sippi_kernels_%s',txt{i}));
    
    
    figure_focus(1)
    subplot(1,7,i+1);
    imagesc(prior{1}.x,prior{1}.y,K{i});axis image
    title(txt{i})
    caxis([-1 1].*.02)
end
colormap(jet)
colorbar_shift;
print_mul('sippi_kernels');
