% nmo_compare_forward
%
% see also nmo_setup_example
%
clear all;

% load data from nmo_setup_example
load nmo_reference_data_type1_SN1
%load nmo_reference_data_type2_SN1


%% TEST DIFFERENT FORWARD MODELS
clear d f type
it=0;
it=it+1;type{it}='akirichards';
it=it+1;type{it}='shuey_2_term';
it=it+1;type{it}='shuey_3_term';
%it=it+1;type{it}='shuey';
it=it+1;type{it}='shuey_castagna';
it=it+1;type{it}='zoeppritz';
it=it+1;type{it}='weak_contrast';

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
  subplot(2,4,itype)
  hist(d{itype}{1}-data_ref{1}.d_obs)
  xlabel('difference')
  title(forward.type);
  
end
