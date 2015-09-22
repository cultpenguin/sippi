% nmo_compare_forward
%
% see also nmo_setup_example
%
clear all;close all;

% load data from nmo_setup_example
mat_file='nmo_reference_data_type2_SN1_nx11.mat';
if exist(mat_file);
  load(mat_file)
else
  SN=1;
  ptype=2;
  nx=11;
  nmo_setup_example
end
%% TEST DIFFERENT FORWARD MODELS
clear d f type
it=0;
it=it+1;type{it}='zoeppritz';
it=it+1;type{it}='akirichards';
it=it+1;type{it}='shuey_2_term';
it=it+1;type{it}='shuey_3_term';
it=it+1;type{it}='shuey_castagna';
it=it+1;type{it}='buland_omre';

data_ref=data;
for i=1:length(data)
  data_ref{i}.d_obs=data_ref{i}.d_ref;
end


for itype=1:length(type)
  forward.type=type{itype};
  [d{itype},f{itype}]=sippi_forward(m_ref,forward,prior);
  
  figure(1);
  subplot(2,4,itype)
  hold on
  sippi_plot_data_reflection_nmo(d{itype},data,1,prior)
  title(forward.type);
  
  figure(2);
  subplot(2,4,itype)
  hold on
  sippi_plot_data_reflection_nmo(d{itype},data_ref,1,prior)
  title(forward.type);
  
  figure(3);
  xlim=[-.05 .05];
  subplot(2,4,itype)
  hist(d{itype}{1}-data_ref{1}.d_obs,linspace(xlim(1),xlim(2),31))
  set(gca,'xlim',xlim);
  xlabel('difference')
  title(forward.type);
  
end

%% TEST MODELING ERROR
prior{1}.d_target=[7.9 8.1];prior{1}.m0=0;
prior{2}.d_target=[7.9 8.1];prior{2}.m0=0;

forward_1=forward;
forward_1.type='zoeppritz';
forward_2=forward;
forward_2.type='weak_contrast';
%forward_2.type='akirichards';
[Ct,dt,dd]=sippi_compute_modelization_forward_error(forward_1, forward_2, prior,100);
%
figure(5);
subplot(2,1,1);
wiggle(forward.angle,forward.t,reshape(dd{1}(:,1),201,10))
title('Modeing error realization')
%subplot(2,2,2);
n=diag(diag(Ct{1})).*.00001;
r=gaussian_simulation_cholesky(dt{1},Ct{1}+n,1);
hold on
wiggle(forward.angle,forward.t,reshape(r,201,10))
hold off

subplot(2,2,3);
plot(sqrt(diag(Ct{1})))
ylabel('std(diag(Ct))')

subplot(2,2,4);
imagesc(Ct{1})
title('Ct_{est}')
caxis([-1 1].*1e-4)

SN=1./(std(reshape(dd{1}(:,1),201,10))./std(reshape(d{1}{1},201,10)))

