% sippi_AM13 Compare traveltime forward models in SIPPI
%
clear all;close all

%% Mke some choices
if ~exist('use_prior','var')
    % use_prior=1; % Gaussian
    % use_prior=2; % Gaussian with bimodal target distribution
    % use_prior=3; % Gaussian with uniform target distribution
    use_prior=4; % Plurigaussian
    % use_prior=5; % Gaussian with variable covariance parameters
    % use_prior=6; % Matern type covariance with varying nu parameter
end
if ~exist('use_forward','var')
    use_forward=1; % ray_2d - linear straight ray
    % use_forward=2; % linear ray (using eikonal)
    % use_forward=3; % bended ray
    % use_forward=4; % linear fat
    % use_forward=5; % bended fat
    % use_forward=6; % eikonal
    % use_forward=7; % waveform FD + first arriavle
end

if ~exist('dx','var')
    dx=0.2;
end


%% SETUP DATA, PRIOR and FORWARD

%% SETUP DATA
D=load('AM13_data.mat');

id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std.*0+0.4;;
%data{id}.i_use=1:20;
%data{id}.Ct=1; % Data covariance describing modelization error
data{id}.Ct=D.Ct; % Correlated noise model according to Cordua et al (2008; 2009)

%% SETUP DIFFERENT PRIOR STRUCTURES
% define some standard values
m0=0.145;
Va='.0003 Sph(6,90,.3)';

% some parameters needed by all a priori types
np=2;
prior_ref{1}.name='Velocity (m/ns)';
prior_ref{1}.x=[(0-np*dx):dx:(5+np*dx)];
prior_ref{1}.y=[(1-np*dx):dx:(12+np*dx)];
prior_ref{1}.cax=[.10 .17];

% define a number of a priori models
im_all=0;

%% GAUSSIAN
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Gaussian';
prior_all{im_all}{im}.type='FFTMA';
prior_all{im_all}{im}.m0=m0;
prior_all{im_all}{im}.Va=Va;

%% GAUSSIAN WITH BIMODAL TARGET DISTRIBUTION
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Gaussian Bimodal';
prior_all{im_all}{im}.type='FFTMA';
prior_all{im_all}{im}.Cm=Va;
% bimodal distribution
N=10000;
prob_chan=0.5;
dd=.015;
d1=randn(1,ceil(N*(1-prob_chan)))*.0025+0.145-dd;  %0.1125;
d2=randn(1,ceil(N*(prob_chan)))*.0025+0.145+dd; %0.155;
d_target=[d1(:);d2(:)];
prior_all{im_all}{im}.d_target=d_target;

%% GAUSSIAN WITH UNIFORM TARGET DISTRIBUTION
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Gaussian Uniform';
prior_all{im_all}{im}.type='FFTMA';
prior_all{im_all}{im}.Cm=Va;
N=10000;
prior_all{im_all}{im}.d_target=(rand(1,N)-.5)*.09+m0;

%% PLURIGAUSSIAN
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Plurigaussian';
prior_all{im_all}{im}.type='plurigaussian';
% 1D plurigaussian (truncated)
prior_all{im_all}{im}.pg_prior{1}.Cm=' 1 Gau(5,90,.5)';
prior_all{im_all}{im}.pg_map=[0.11 .11 .13 .13 .15 .15 .17 .13 .17];
% 2D plurigaussian
prior_all{im_all}{im}.pg_prior{2}.Cm=' 1 Sph(10,45,.2)';
prior_all{im_all}{im}.pg_map=[0.11 .11 .17 ; .13 .15 .17 ; .17 .13 .17];

%% GAUSSIAN with variable covariance model parameters
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Gaussian - variable pars';

prior_all{im_all}{im}.type='fftma';
prior_all{im_all}{im}.Va=Va;
prior_all{im_all}{im}.m0=m0;
i_master=im;

% range - horizontal
im=im+1;
prior_all{im_all}{im}.type='uniform';
prior_all{im_all}{im}.name='range_1';
prior_all{im_all}{im}.min=1;
prior_all{im_all}{im}.max=6;
prior_all{im_all}{im}.prior_master=i_master;

% range - vertical
im=im+1;
prior_all{im_all}{im}=prior_all{im_all}{im-1};
prior_all{im_all}{im}.name='range_2';

% rotation
im=im+1;
prior_all{im_all}{im}.type='gaussian';
prior_all{im_all}{im}.name='ang_1';
prior_all{im_all}{im}.m0=90;
prior_all{im_all}{im}.std=10;
prior_all{im_all}{im}.norm=2;
prior_all{im_all}{im}.prior_master=i_master;

% m0
im=im+1;
prior_all{im_all}{im}.type='uniform';
prior_all{im_all}{im}.name='m0';
prior_all{im_all}{im}.min=0.10;
prior_all{im_all}{im}.max=0.17;
prior_all{im_all}{im}.prior_master=i_master;

% sill
im=im+1;
prior_all{im_all}{im}.type='uniform';
prior_all{im_all}{im}.name='sill';
prior_all{im_all}{im}.min=0.001;
prior_all{im_all}{im}.max=0.006;
prior_all{im_all}{im}.prior_master=i_master;

%% MATERN TYPE
im_all=im_all+1;
im=1;
prior_all{im_all}{im}=prior_ref{1};
prior_all{im_all}{im}.title='Gaussian - Matern';

prior_all{im_all}{im}.type='fftma';
prior_all{im_all}{im}.m0=m0;
prior_all{im_all}{im}.Cm='.0003 Mat(2,90,.3)';
prior_all{im_all}{im}.fftma_options.pad_x=150; % avoid striping
prior_all{im_all}{im}.fftma_options.pad_y=150; % avoid striping
i_master=im;

% NU / Matern - variable parameter
im=im+1;
prior_all{im_all}{im}.type='uniform';
prior_all{im_all}{im}.name='nu';
prior_all{im_all}{im}.min=.1;
prior_all{im_all}{im}.max=2;
prior_all{im_all}{im}.prior_master=i_master;

%% PLOT SAMPLE FROM PRIOR_MUL
do_plot_prior_mul=0;
if do_plot_prior_mul==1;
    figure(2);clf;
    nsim=6;
    np=length(prior_all);
    for i=1:np
        prior=prior_all{i};
        
        for isim=1:nsim
            
            subplot(np,nsim,(i-1)*nsim+isim);
            [m,prior]=sippi_prior(prior);
            imagesc(prior{1}.x, prior{1}.y, m{1});
            set(gca,'FontSize',6)
            axis image;
            caxis(prior{1}.cax);
            if isim==1;
                t=title(sprintf('%s) %s',char(i+96),prior{1}.title));
                pos=get(t,'position');
                pos(1)=-5;
                set(t,'position',pos,'HorizontalAlignment','Left','FontSize',10)
            end
        end
        colorbar_shift;
        
    end
    print_mul(sprintf('%s_prior_reals',mfilename))
end


%% SETUP THE FORWARD MODEL(S)
forward_ref.forward_function='sippi_forward_traveltime';
forward_ref.sources=D.S;
forward_ref.receivers=D.R;



k=0;
k=k+1; % RAY
f{k}.type='ray_2d';
f{k}.r=1;
f{k}.name=f{1}.type;

k=k+1; 
f{k}.type='ray';
f{k}.linear=1;
f{k}.name='SR';

k=k+1; 
f{k}.type='ray';
f{k}.linear=0;
f{k}.name='BR';

k=k+1; 
f{k}.type='fat';
f{k}.linear=1;
f{k}.freq=0.1;
f{k}.name='SF';

k=k+1; 
f{k}.type='fat';
f{k}.linear=0;
f{k}.freq=0.1;
f{k}.name='BF';

k=k+1;
f{k}.type='eikonal';
f{k}.name=f{k}.type;

%k=k+1;
%f{k}.type='fd';
%f{k}.name=f{k}.type;

%%

% select prior type
prior=prior_all{use_prior};
% generate a realization from the prior
[m,prior]=sippi_prior(prior);

for k=1:length(f);
    % set forward model
    forward = forward_ref;
    fn=fieldnames(f{k})
    for i=1:length(fn); forward.(fn{i}) = f{k}.(fn{i});end
    disp(sprintf('working on %s',f{k}.name))
    
    [d,forward]=sippi_forward(m,forward,prior,data);
    t_all(:,k)=d{1};
    L{k}=f{k}.name;
    
end

%% SAVE DATA
M = [forward.sources,forward.receivers,t_all]
save('AM13_SR_traveltime.dat','-ascii','M')

V = m{1};
save('AM13_V.dat','-ascii','V')

%% PLOT DATA
figure(1);
subplot(1,2,1);
imagesc(prior{1}.x,prior{1}.y,m{1});
axis image;colorbar
xlabel('X (m)'); ylabel('Y (m)')
subplot(1,2,2);

plot(t_all);
legend(L)
xlabel('Data #')
ylabel('Traveltime (ns)')
print_mul(sprintf('AM13_compare forward'))




