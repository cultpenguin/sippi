% sippi_forward_traveltime_example: Different examples of travel time
%      computation
%
% See also: sippi_forward_traveltime
%
%% load some data
clear all;close all
rng(1)
D=load('AM13_data.mat');
% SETUP DATA
id=1;
data{id}.d_obs=D.d_obs;
data{id}.d_std=D.d_std;


%% DEFINE PRIOR
im=1;
prior{im}.type='FFTMA';
prior{im}.name='Velocity (m/ns)';
prior{im}.Va='.0003 Sph(6,90,.33)';
dx=0.075;
prior{im}.x=[-1:dx:6];
prior{im}.y=[0:dx:13];
prior{im}.cax=[.04 .18];

d_target=[randn(1,100)*.003+0.11 randn(1,100)*.003+0.16];
prior{im}.d_target=d_target;


%% SETUP THE FORWARD MODEL
forward.forward_function='sippi_forward_traveltime';
forward.sources=D.S;
forward.receivers=D.R;
[m,prior]=sippi_prior(prior);
sippi_plot_prior(prior,m);

for it=1:4;
    forwardi{it}=forward;
    if it==1;
        forwardi{it}.type='fd';
        %forwardi{it}.fa.doPlot=0;
        %forwardi{it}.fa.use_method=2;
        forwardi{it}.freq=0.1;
    elseif it==2
        forwardi{it}.type='eikonal';
    elseif it==3
        forwardi{it}.type='fat';
        forwardi{it}.linear=1;
        forwardi{it}.freq=0.1;
    elseif it==4
        forwardi{it}.type='fat';
        forwardi{it}.linear=0;
        forwardi{it}.freq=0.1;
        %forwardi{it}.linear_m=m{1};
    end

    L{it}=forwardi{it}.type;
    [d{it},forwardi{it},prior,data]=sippi_forward_traveltime(m,forwardi{it},prior,data);
    
    figure(3);
    plot(d{it}{1});
    hold on
end
hold off
legend(L)