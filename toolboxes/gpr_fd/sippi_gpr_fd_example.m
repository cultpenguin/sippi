% sippi_setup_gpr_fd: Inversion of full waveform GPR data.
%
%
% /REF/See also
%    https://sites.google.com/site/kscordua/downloads/tomographic-full-waveform-inversion
%    Ernst, J. R., H. Maurer, A. G. Green, and K. Holliger, 2007, Full-waveform inversion of crosshole radar data based on 2-D finite-difference time domain solutions of Maxwell�s equations: IEEE Transactions on Geoscience and Remote Sensing, 45, 2807�2828.
%    Cordua, K. S., T. M. Hansen, and K. Mosegaard, 2012. Monte Carlo full waveform inversion of crosshole GPR data using multiple-point geostatistical a priori information. Geophysics, 77(2), H19 � H31.
%
%% PRIOR 
rng(1)
clear all;
% EPS
ip=1;
dx=0.1;nx=52;ny=120;
dx=0.05;nx=104;ny=240;
prior{ip}.type='fftma';
prior{ip}.name='eps';
prior{ip}.x=[1:nx].*dx;
prior{ip}.y=[1:ny].*dx;

use_target=0;
d_target=[ones(1,1801).*2.5,ones(1,4432).*4.5];
if use_target==0;
  prior{ip}.m0=mean(d_target);
  v0=var(d_target);
  prior{ip}.Cm=sprintf('%4.4g Sph(5,90,.15)',v0);
else
  prior{ip}.Cm='1 Sph(5,90,.15)';
  prior{ip}.d_target=d_target;  
 
end

% SIG

%% FORWARD
ant_pos=[
    0.1 3 5.1 0.1
    0.1 3 5.1 1.5
    0.1 3 5.1 3
    0.1 3 5.1 4.5
    0.1 3 5.1 6
    0.1 9 5.1 6
    0.1 9 5.1 7.5
    0.1 9 5.1 9
    0.1 9 5.1 10.5
    0.1 9 5.1 11.9
    5.1 3 0.1 0.1
    5.1 3 0.1 1.5
    5.1 3 0.1 3
    5.1 3 0.1 4.5
    5.1 3 0.1 6
    5.1 9 0.1 6
    5.1 9 0.1 7.5
    5.1 9 0.1 9
    5.1 9 0.1 10.5
    5.1 9 0.1 11.9];

forward.ant_pos=ant_pos;
forward.sig=3; % constant
forward.output_type='trace'; % 'trace' or 'gather' 
forward.output_it=10; % output every 'output_ti' samples
forward.t=1e-7; % SIMULATION TIME
forward.sources=ant_pos(:,1:2);
forward.receivers=ant_pos(:,3:4);
    

% options for forward modelerreturn
forward.addpar.debug=-1;
forward.addpar.cores=4;
forward.addpar.Tg=100*10^6;
%forward.addpar.executable='/path/to/executable';
%forward.addpar.work_dir='/path/to/working_directory';


[m_ref,prior]=sippi_prior(prior);
sippi_plot_prior(prior,m_ref);
[d,forward,prior]=sippi_forward_gpr_fd(m_ref,forward,prior);
sippi_plot_data_gpr(d);
return


%% TEST FORWARD / MAKE REFERENCE DATA
[m_ref,prior]=sippi_prior(prior);
sippi_plot_prior(prior,m_ref);
[d,forward,prior]=sippi_forward_gpr_fd(m_ref,forward,prior);



return
forward_t=forward;
forward_t.type='eikonal';
[d_t,forward_t,prior]=sippi_forward_traveltime(m_ref,forward_t,prior);


for id=1:length(d);
  data{id}.d_ref=d{id};
  data{id}.d_std=.001;
  
  data{id}.d_noise=randn(size(data{id}.d_ref)).*data{id}.d_std;
  data{id}.d_obs=data{id}.d_ref+data{id}.d_noise;
  
end
sippi_plot_data_gpr(d,data);
print -dpng ref_data;
figure(10);imagesc(m_ref{1});axis image;print -dpng ref_model
return
%% Metropolis sampling
%forward.addpar.debug=-1; % NO OUTPUT
forward.forward_function='sippi_forward_gpr_fd';
prior{1}.seq_gibbs.step=1;
prior{1}.seq_gibbs.i_update_step=25;
options.mcmc.m_ref=m_ref;
options.mcmc.nite=20000;   % [1] : Number if iterations
options.mcmc.i_sample=50; % : Number of iterations between saving model to disk
options.mcmc.i_plot=50;  % [1]: Number of iterations between updating plots
options.mcmc.i_save_workspace=1000;  % [1]: Number of iterations between complete save of workspace
[options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options);
save(options.txt)

sippi_plot_posterior(options.txt);


