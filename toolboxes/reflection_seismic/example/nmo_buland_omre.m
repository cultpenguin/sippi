% nmo_compare_forward
%
% see also nmo_setup_example
%
clear all;close all;

% load data from nmo_setup_example
mat_file='nmo_reference_data_type2_SN10_nx11.mat';
if exist(mat_file);
  load(mat_file)
else
  SN=10;
  ptype=2;
  nx=11;
  nmo_setup_example
end

%% Perform linearized inversion

forward.type='weak_contrast';
[d,f]=sippi_forward(m_ref,forward,prior);
%setup Cm
for i=[1,3,2];
  pos=[x(1)+t(:).*0 t(:)];
  C=precal_cov(pos,pos,prior{i}.Va);
  if i==1;
    Cm=C;
  else
    Cm=blkdiag(Cm,C);
  end
end
%Cm=eye(3*nt).*var_omre(1);
  
 [M,V,vp_est,vs_est,rho_est]=buland_omre_inversion(d,f.vp0,f.vs0,f.rho0,f.angle,f.wl,Cd,Cm);

 %%
 figure(9);
 subplot(2,4,1);
 imagesc(m_ref{1});cax=caxis;
 title('Vp true')
 subplot(2,4,5);
 imagesc(vp_est);caxis(cax);
 title('Vp lsq est')
 
 subplot(2,4,2);
 imagesc(m_ref{2});cax=caxis;
 title('Vs')
 subplot(2,4,6);
 imagesc(vs_est);caxis(cax);
 
 subplot(2,4,3);
 imagesc(m_ref{3});cax=caxis;
 title('rho')
 subplot(2,4,7);
 imagesc(rho_est);caxis(cax);
 
 subplot(2,4,4);
 imagesc(m_ref{3}.*m_ref{1});cax=caxis;
 title('IP')
 subplot(2,4,8);
 imagesc(rho_est.*vp_est);caxis(cax);
 
 %% 
 figure(10);
 ix=1;
 t=forward.t;
 nt=length(t);
 
 subplot(1,4,1);
 plot(exp(vp_est(:,ix)),t,'k-')
 d_std=sqrt(V([1:nt],ix));
 hold on
 plot(exp(vp_est(:,ix)+2*d_std),t,'k--')
 plot(exp(vp_est(:,ix)-2*d_std),t,'k--')
 plot(exp(m_ref{1}(:,ix)),t,'r-')
 hold off
 title('Vp')
 
 subplot(1,4,2);
 plot(exp(vs_est(:,ix)),t,'k-')
 d_std=sqrt(V([1:nt]+nt,ix));
 hold on
 plot(exp(vs_est(:,ix)+2*d_std),t,'k--')
 plot(exp(vs_est(:,ix)-2*d_std),t,'k--')
 plot(exp(m_ref{2}(:,ix)),t,'r-')
 hold off
 title('Vs')
 
 subplot(1,4,3);
 plot(exp(rho_est(:,ix)),t,'k-')
 d_std=sqrt(V([1:nt]+nt,ix));
 hold on
 plot(exp(rho_est(:,ix)+2*d_std),t,'k--')
 plot(exp(rho_est(:,ix)-2*d_std),t,'k--')
 plot(exp(m_ref{3}(:,ix)),t,'r-')
 hold off
 title('rho')
 
 
 figure(12);
 vp=M(1:201,:);
 plot(m_ref{1},vp,'.')
 
 
  