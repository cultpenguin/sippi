clear all,close all

for ij=1:5
    
    if ij==1
        file='wedge_reference_data_corr_0200_m3.mat';
    elseif ij==2
        file='wedge_reference_data_corr_0100_m3.mat';
    elseif ij==3
        file='wedge_reference_data_corr_0050_m3.mat';
    elseif ij==4
        file='wedge_reference_data_corr_0020_m3.mat';
    elseif ij==5
        file='wedge_reference_data_corr_0010_m3.mat';
    end
    fp='..\..\..\..\SIPPI_CERE\examples\wedge\REFERENCE_DATA';
    fp='DATA';
    
    MAT=load([fp,filesep,file]);
    t=MAT.forward.t;
    nt=MAT.forward.nt;
    nx=MAT.forward.nx;
    d_ref=MAT.data{1}.d_obs;
    
    figure(1);
    imagesc(reshape(d_ref,nt,nx))
    iang=1;
    Cd=MAT.data{iang}.Cd;
    %Cd=diag(diag(Cd)); % IGNORE CORRELATED NOISE
    wl=MAT.forward.wl(iang);ang=MAT.forward.angle(iang);
    vp0=2000+250;
    %vp0=2000; % BASE
    vs0=vp0/sqrt(2);
    rho0=vp0;
    Cm_old=load('Cm');;
    Cm_unit=eye(nt,nt);
    useCmEye=0;
    if useCmEye==1
        Va=sprintf('1 Sph(%f)',.01*MAT.forward.dt);
    else
        Va=sprintf('1 Sph(%f)',30*MAT.forward.dt);
    end
    Cm_unit=precal_cov(MAT.forward.t(:),MAT.forward.t(:),Va);
    var_arr=[0.0074 0.0024 0.0074];
    var_arr=[1 1/1.4 1].*0.0124;
    Cm=blkdiag(Cm_unit.*var_arr(1),Cm_unit.*var_arr(2),Cm_unit.*var_arr(3));
    
    
    %% SETUP DATA
    for id=1:nx;
        it=[1:nt]+(id-1)*nt;
        d{id}=d_ref(it);
    end
    
    [M,V,vp_est,vs_est,rho_est,Cm_est]=buland_omre_inversion(d,log(vp0),log(vs0),log(rho0),ang,wl,Cd,Cm);
    
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
    
    %% VP PROB
    M_VP=M(1:nt,:);
    V_VP=V(1:nt,:);
    V_thres=vp0+200;
    P = 1-normcdf(log(V_thres),M_VP,sqrt(V_VP));
    figure(5);clf
    imagesc(0:40,t,P(:,20:60))
    xlabel('X')
    ylabel('Twi way time (s)')
    colormap(hot)
    ppp(10,6,4,2,2)
    set(gca,'FontSize',6)
    caxis([0 1])
    print_mul(sprintf('prob_layer_lsq_Cm%d_%s',useCmEye,file(22:30)))
    
    %% PRIOR REALIZATIONS
    nr=51;
    vp_prior=gaussian_simulation_cholesky(log(vp0),Cm(1:nt,1:nt),31);
    d50{1}=d{60};
    [Ms,Vs,vp_ests,vs_ests,rho_ests,Cm_est]=buland_omre_inversion(d50,log(vp0),log(vs0),log(rho0),ang,wl,Cd,Cm);
    vp_post=gaussian_simulation_cholesky(Ms(1:nt),Cm_est(1:nt,1:nt),31);
    figure(11);set_paper('portrait');clf;
    plot(exp(vp_prior),t,'k-');
    xlabel('Vp (m/s)')
    ylabel('Two way time (s)')
    axis([1500 3300 0.15 0.30])
    set(gca,'ydir','reverse')
    ppp(5,8,10,3,2)
    print_mul(sprintf('prob_layer_reals_prior_Cm%d_%s',useCmEye,file(22:30)))
    hold on;
    plot(exp(vp_post),t,'r-');
    plot([2000 2000 2500 2500 2000 2000],[0 0.2 0.2 0.25 0.25 0.4],'g--','LineWidth',3)
    hold off
    print_mul(sprintf('prob_layer_reals_post_Cm%d_%s',useCmEye,file(22:30)))
end
%%  PRIOR PROB
M_VP=M(1:nt,:).*0+log(vp0);
V_VP=V(1:nt,:).*0+var_arr(1);;
V_thres=vp0+200;
P = 1-normcdf(log(V_thres),M_VP,sqrt(V_VP));
figure(6);clf
imagesc(0:40,t,P(:,20:60))
xlabel('X')
ylabel('Two way time (s)')
colormap(hot)
ppp(10,6,4,2,2)
set(gca,'FontSize',6)
caxis([0 1])
print_mul(['prob_layer_lsq_prior'])
%colorbar
%print_mul(['prob_layer_lsq_prior_cbar'])


%% REALS FROM OPTIMAL PRIOR


figure(11);set_paper('portrait');clf;
for ir=1:nr,;
    t1=ceil(rand(1)*401);
    t2=t1+ceil(rand(1)*(401-t1));
    plot(ir+[2000 2000 2500 2500 2000 2000],[0 t1./1000 t1./1000 t2./1000 t2./1000 0.4],'k-')
    hold on
end
xlabel('Vp (m/s)')
ylabel('Two way time (s)')
axis([1500 3300 0.15 0.30])
set(gca,'ydir','reverse')
ppp(5,8,10,3,2)
print_mul(sprintf('prob_layer_reals_prior_Cm%d_optimal',useCmEye))
hold on;
plot([2000 2000 2500 2500 2000 2000],[0 0.2 0.2 0.25 0.25 0.4],'r-','LineWidth',4)
plot([2000 2000 2500 2500 2000 2000],[0 0.2 0.2 0.25 0.25 0.4],'g--','LineWidth',3)
hold off
print_mul(sprintf('prob_layer_reals_post_Cm%d_optimal',useCmEye))





