clear all,close all

file='wedge_reference_data_corr_0400_m3.mat';
file='wedge_reference_data_corr_0100_m3.mat';
file='wedge_reference_data_corr_0050_m3.mat';
%file='wedge_reference_data_corr_0010_m3.mat';
fp='..\..\..\..\SIPPI_CERE\examples\wedge\REFERENCE_DATA';

MAT=load([fp,filesep,file]);

nt=MAT.forward.nt;
nx=MAT.forward.nx;
d_ref=MAT.data{1}.d_obs;

figure(1);
imagesc(reshape(d_ref,nt,nx))
iang=1;
Cd=MAT.data{iang}.Cd;
wl=MAT.forward.wl(iang);ang=MAT.forward.angle(iang);
%vp0=2000+500;
vp0=2000; % BASE
vs0=vp0/sqrt(2);
rho0=vp0;
Cm_old=load('Cm');;
Cm_unit=eye(nt,nt);
%Va=sprintf('1 Sph(%f)',.01*MAT.forward.dt);
Va=sprintf('1 Sph(%f)',30*MAT.forward.dt);
Cm_unit=precal_cov(MAT.forward.t(:),MAT.forward.t(:),Va);
var_arr=[0.0074 0.0024 0.0074];
Cm=blkdiag(Cm_unit.*var_arr(1),Cm_unit.*var_arr(2),Cm_unit.*var_arr(3));


%% SETUP DATA
for id=1:nx;
    it=[1:nt]+(id-1)*nt;
    d{id}=d_ref(it);
end

[M,V,vp_est,vs_est,rho_est]=buland_omre_inversion(d,log(vp0),log(vs0),log(rho0),ang,wl,Cd,Cm);

%%
ip_est=exp(vp_est).*exp(rho_est);
figure(2);
imagesc(MAT.forward.x,MAT.forward.t,ip_est);

figure(3);
subplot(1,3,1);
plot(ip_est(:,60))
subplot(1,3,2);
plot(exp(vp_est(:,60)))
subplot(1,3,3);
plot(exp(rho_est(:,60)))
