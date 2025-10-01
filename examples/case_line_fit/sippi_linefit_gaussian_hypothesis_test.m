% sippi_linefit: Fiting line using SIPPI
clear all;
%rng('default')
rng(randi(10000))

%% Load or generate data data
nd=11;
%
mat_datfile = sprintf('sippi_linefit_data_%d',nd);
if exist([mat_datfile,'.mat'])
    load(mat_datfile);
else
    sippi_linefit_make_data;
end

x_obs = x;

x = 0:1:30;

%% setup the data
%d_std = 2*d_std;
d_var = d_std.^2;
data{1}.d_obs=d_obs;
data{1}.d_std=d_std;

%% setup the forward
forward.x=x;
forward.forward_function='sippi_forward_linefit';

for ix=1:length(x_obs)
    ix_obs(ix)=find(x_obs(ix)==forward.x);
end

%% setup the prior model
% the intercept
im=1;
prior{im}.type='gaussian';
prior{im}.name='intercept';
prior{im}.m0=0;
prior{im}.std=30;
prior{im}.m_true=m_ref{1};

% 1st order, the gradient
im=2;
prior{im}.type='gaussian';
prior{im}.name='gradient';
prior{im}.m0=0;
prior{im}.std=4;
prior{im}.norm=80;
prior{im}.m_true=m_ref{2};

% 2nd order
im=3;
prior{im}.type='uniform';
prior{im}.name='2nd';
prior{im}.min=-.1;
prior{im}.max=.1;
prior{im}.m_true=m_ref{3};

% 3rd order
im=4;
prior{im}.type='uniform';
prior{im}.name='4rd';
prior{im}.min=-.003;
prior{im}.max=.003;

p{1}{1}=prior{1};

p{2}{1}=prior{1};
p{2}{2}=prior{2};

p{3}{1}=prior{1};
p{3}{2}=prior{2};
p{3}{3}=prior{3};

p{4}{1}=prior{1};
p{4}{2}=prior{2};
p{4}{3}=prior{3};
p{4}{4}=prior{4};

%% TEST
clear d_sim*
N=1000000;
np=length(p);
for ip=1:np
    progress_txt(ip,np)
    for i=1:N
        m=sippi_prior(p{ip});
        d=sippi_forward(m,forward);

        d_sim{ip}(i,:)=d{1};
        d_sim_obs{ip}(i,:)=d{1}(ix_obs);
        dd = d_sim_obs{ip}(i,:)'-d_obs;
        logL_sim{ip}(i)= -0.5*sum(dd.^2./d_var);

        %plot(forward.x,d{1},'-',Color=col)
        %hold on
    end
    EV(ip)=log(mean(exp(logL_sim{ip})));
    logL_max(ip) = max(logL_sim{ip});
    P_acc{ip}=exp(logL_sim{ip}-logL_max(ip));
    r_acc = rand(1,N);
    i_acc{ip}=find(P_acc{ip}>r_acc);
    n_acc = length(i_acc{ip});

end

P_hyp=exp(EV)/sum(exp(EV));

%% INVERSION ON MIXTURE
n_obs=length(d_obs);
d_sim_obs_mix = zeros(np*N,n_obs);
hypothesis_mix= zeros(np*N,1);
for ip=1:np
    i1=(ip-1)*N+1;
    i2=i1-1+N;
    d_sim_obs_mix(i1:i2,:)=d_sim_obs{ip};
    d_sim_mix(i1:i2,:)=d_sim{ip};
    hypothesis_mix(i1:i2)=ip;
end
for i=1:(np*N)
    if mod(i,N)==0;progress_txt(i,np*N);end
    dd = d_sim_obs_mix(i,:)'-d_obs;
    logL_sim_mix(i)= -0.5*sum(dd.^2./d_var);
end
EV_mix=log(mean(exp(logL_sim_mix)));
logL_max_mix = max(logL_sim_mix);
P_acc_mix=exp(logL_sim_mix-logL_max_mix);
r_acc_mix = rand(1,np*N);
i_acc_mix=find(P_acc_mix>r_acc_mix);
n_acc = length(i_acc_mix);

for ip=1:length(p)
    n_c(ip)=length(find(hypothesis_mix(i_acc_mix)==ip));
end
P_hyp_mix=n_c/sum(n_c);

disp(['P_hyp (E)    = [',sprintf('%3.3f ',P_hyp),']'])
disp(['P_hyp_mix = [',sprintf('%3.3f ',P_hyp_mix),']'])


%%
figure(1);clf
bar([P_hyp;P_hyp_mix]')
legend({'Evidence','Mixed prior'})
xlabel('Prior type')
ylabel('P(prior|data)')
ylim([0,.5])
print_mul(sprintf('EV_compare_N%d_%s',N,datestr(now,'HH_MM_SS')))
%%

figure(2);clf
for ip=1:length(p);
    col= [0 .9 .5].*ip/np;
    col= [0 0 0];
    nshow=min([10,length(i_acc{ip})]);
    plot(forward.x,d_sim{ip}(i_acc{ip}(1:nshow),:),'-',Color=col)
    hold on
end
nshow=min([40,length(i_acc_mix)]); 
i_acc_mix=shuffle(i_acc_mix);
plot(forward.x,d_sim_mix(i_acc_mix(1:nshow),:),'b-')
errorbar(x_obs,d_obs,2*d_std,'k*')
hold off
print_mul(sprintf('forward_N%d_%s',N,datestr(now,'HH_MM_SS')))

%% SAVE TO MAT FILES
doSave=true;
if doSave
    disp('saving output')
    save(sprintf('linefit_evidence_N%d_std%d_%s',N,d_std(1),datestr(now,'HH_MM_SS')),'P_hyp', ...
        'P_hyp_mix', ...
        'd_sim_obs_mix', 'd_sim_obs', ...
        'd_sim_mix', 'd_sim', ...
        'd_obs','d_std', ...
        'logL_sim_mix','logL_sim', ...
        'd_sim_obs','d_sim_obs_mix', 'hypothesis_mix','EV','EV_mix')
end