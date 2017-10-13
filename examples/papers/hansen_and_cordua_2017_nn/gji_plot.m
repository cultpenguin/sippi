loadData=1;
if loadData==1;
    
    clear all;close all;
    
    try   
        load('grl_ReferenceModel','m_ref');
        load('grl_ReferenceModel','prior');
        load('grl_ReferenceModel','forward');
        load('grl_ReferenceModel','data');
    catch
        load('grl_ReferenceModel_t0','m_ref');
        load('grl_ReferenceModel_t0','prior');
        load('grl_ReferenceModel_t0','forward');
        load('grl_ReferenceModel_t0','data');
    end
    %load grl_NM2376_DX20_fd_NT40000_SD3_NH80_modelerr;
    %load grl_NM2376_DX20_fd_NT40000_SD3_NH80_inverted;
    txt='grl';
    
    d_txt=dir('grl_*_inverted.mat');
    disp(sprintf('loading data: %s',d_txt(1).name));
    load(d_txt(1).name);
    
    clear tit o 
    % SETUP tit and o
    if exist('f_mul');
        for io=1:length(f_mul);
            o{io}=options_out{io}.txt;
            if strcmp(f_mul{io}.forward_function,'sippi_forward_mynn')
                tit{io}=sprintf('g_{%d}',size(f_mul{io}.ATTS,2));                
            elseif strcmp(f_mul{io}.forward_function,'sippi_forward_traveltime')
                tit{io}=sprintf('g_{%s}',f_mul{io}.type(1:3));
            else
                tit{io}=sprintf('g_{%d}',io);
            end
        end
    end
    
    
    
    %i=i+1;o{i}='20161211_2025_sippi_metropolis_grl_NM2376_DX20_fd_NT1000_SD3_NH80';Nf_arr(i)=1000;tit{i}='1000';
    %i=i+1;o{i}='20161211_2037_sippi_metropolis_grl_NM2376_DX20_fd_NT5000_SD3_NH80';Nf_arr(i)=5000;;tit{i}='5000';
    %i=i+1;o{i}='20161211_2050_sippi_metropolis_grl_NM2376_DX20_fd_NT10000_SD3_NH80';Nf_arr(i)=10000;tit{i}='10000';
    %i=i+1;o{i}='20161211_2103_sippi_metrop
    
    %i=i+1;o{i}='20161211_2025_sippi_metropolis_grl_NM2376_DX20_fd_NT1000_SD3_NH80';Nf_arr(i)=1000;tit{i}='1000';
    %i=i+1;o{i}='20161211_2037_sippi_metropolis_grl_NM2376_DX20_fd_NT5000_SD3_NH80';Nf_arr(i)=5000;;tit{i}='5000';
    %i=i+1;o{i}='20161211_2050_sippi_metropolis_grl_NM2376_DX20_fd_NT10000_SD3_NH80';Nf_arr(i)=10000;tit{i}='10000';
    %i=i+1;o{i}='20161211_2103_sippi_metropolis_grl_NM2376_DX20_fd_NT20000_SD3_NH80';Nf_arr(i)=20000;tit{i}='20000';
    %i=i+1;o{i}='20161211_2117_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80';Nf_arr(i)=40000;tit{i}='40000';
    %i=i+1;o{i}='20161211_2133_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80_eikonal';Nf_arr(i)=0;tit{i}='eik';
    %i=i+1;o{i}='20161211_2237_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80_ray_2d';Nf_arr(i)=0;tit{i}='g_{ray}';
    
    % LONG RUN nCt=6000
    %i=i+1;o{i}='20161213_1524_sippi_metropolis_grl_NM2376_DX20_fd_NT1000_SD3_NH80';Nf_arr(i)=1000;tit{i}='g_{1000}';
    %i=i+1;o{i}='20161213_1718_sippi_metropolis_grl_NM2376_DX20_fd_NT5000_SD3_NH80';Nf_arr(i)=5000;;tit{i}='g_{5000}';
    %i=i+1;o{i}='20161213_1917_sippi_metropolis_grl_NM2376_DX20_fd_NT10000_SD3_NH80';Nf_arr(i)=10000;tit{i}='g_{10000}';
    %i=i+1;o{i}='20161213_2157_sippi_metropoolis_grl_NM2376_DX20_fd_NT20000_SD3_NH80';Nf_arr(i)=20000;tit{i}='20000';
    %i=i+1;o{i}='20161211_2117_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80';Nf_arr(i)=40000;tit{i}='40000';
    %i=i+1;o{i}='20161211_2133_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80_eikonal';Nf_arr(i)=0;tit{i}='eik';
    %i=i+1;o{i}='20161211_2237_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80_ray_2d';Nf_arr(i)=0;tit{i}='g_{ray}';
    
    % LONG RUN nCt=6000
    %i=i+1;o{i}='20161213_1524_sippi_metropolis_grl_NM2376_DX20_fd_NT1000_SD3_NH80';Nf_arr(i)=1000;tit{i}='g_{1000}';
    %i=i+1;o{i}='20161213_1718_sippi_metropolis_grl_NM2376_DX20_fd_NT5000_SD3_NH80';Nf_arr(i)=5000;;tit{i}='g_{5000}';
    %i=i+1;o{i}='20161213_1917_sippi_metropolis_grl_NM2376_DX20_fd_NT10000_SD3_NH80';Nf_arr(i)=10000;tit{i}='g_{10000}';
    %i=i+1;o{i}='20161213_2157_sippi_metropolis_grl_NM2376_DX20_fd_NT20000_SD3_NH80';Nf_arr(i)=20000;tit{i}='g_{20000}';
    %i=i+1;o{i}='20161214_0023_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80';Nf_arr(i)=40000;tit{i}='g_{40000}';
    %i=i+1;o{i}='20161214_0315_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80_eikonal';Nf_arr(i)=0;tit{i}='g_{eik}';
    %i=i+1;o{i}='20161214_2128_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80_ray_2d';Nf_arr(i)=0;tit{i}='g_{ray}';
    
    %     %LONGEPOCH
    %     i=0;
    %     i=i+1;o{i}='20161215_2324_sippi_metropolis_grl_NM2376_DX20_fd_NT1000_SD3_NH80';Nf_arr(i)=1000;tit{i}='g_{1000}';
    %     i=i+1;o{i}='20161216_0029_sippi_metropolis_grl_NM2376_DX20_fd_NT5000_SD3_NH80';Nf_arr(i)=1000;tit{i}='g_{5000}';
    %     i=i+1;o{i}='20161216_0135_sippi_metropolis_grl_NM2376_DX20_fd_NT10000_SD3_NH80';Nf_arr(i)=1000;tit{i}='g_{10000}';
    %     i=i+1;o{i}='20161216_0240_sippi_metropolis_grl_NM2376_DX20_fd_NT20000_SD3_NH80';Nf_arr(i)=1000;tit{i}='g_{20000}';
    %     i=i+1;o{i}='20161216_0352_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80';Nf_arr(i)=1000;tit{i}='g_{40000}';
    %     i=i+1;o{i}='20161216_0510_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80_eikonal';Nf_arr(i)=1000;tit{i}='g_{eik}';
    %     i=i+1;o{i}='20161216_1219_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80_ray_2d';Nf_arr(i)=1000;tit{i}='g_{ray}';
    %
    %     d_dir=dir('*sippi_metropolis_grl*');
    %     clear o;
    %     for i=1:length(d_dir)
    %         o{i}= d_dir(i).name;
    %
    %     end
    %
    %
    
    
    
    % LONG RUN
    %     i=i+1;o{i}='20161211_2242_sippi_metropolis_grl_NM2376_DX20_fd_NT1000_SD3_NH80';Nf_arr(i)=1000;tit{i}='g_{1000}';
    %     i=i+1;o{i}='20161212_0259_sippi_metropolis_grl_NM2376_DX20_fd_NT5000_SD3_NH80';Nf_arr(i)=5000;;tit{i}='g_{5000}';
    %     i=i+1;o{i}='20161212_0710_sippi_metropolis_grl_NM2376_DX20_fd_NT10000_SD3_NH80';Nf_arr(i)=10000;tit{i}='g_{10000}';
    %     i=i+1;o{i}='20161212_1135_sippi_metropolis_grl_NM2376_DX20_fd_NT20000_SD3_NH80';Nf_arr(i)=20000;tit{i}='g_{20000}';
    %     i=i+1;o{i}='20161212_1554_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80';Nf_arr(i)=40000;tit{i}='g_{40000}';
    %     i=i+1;o{i}='20161212_2044_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80_eikonal';Nf_arr(i)=0;tit{i}='g_{eik}';
    %    i=i+1;o{i}='20161213_1924_sippi_metropolis_grl_NM2376_DX20_fd_NT40000_SD3_NH80_ray_2d';Nf_arr(i)=0;tit{i}='g_{ray}';
    %
    
    % EX TEST40000_6000_NH40
    %load grl_NM2376_DX20_fd_NT40000_SD3_NH40_inverted.mat
    
    % EX Long Epoch
    % load grl_NM2376_DX20_fd_NT40000_SD3_NH80_inverted.mat
    
    %%TEST40000_6000
    %load grl_NM2376_DX20_fd_NT40000_SD3_NH80_inverted.mat
    
    % TEST 4000 (1000)
    %load TEST40000/grl_NM2376_DX20_fd_NT40000_SD3_NH80_inverted.mat
    
    
    
    % READ DATA
    for io=1:length(o);
        
        
        % [r{i},em{i},ev{i},r_all{i}]=sippi_get_sample(o{i});
        [r{io},em{io},ev{io},r_all{io}]=sippi_get_sample(o{io},1,15,0);
        
        MEAN_POST_VAR(io)=mean(ev{io}(:));
        
        m=load([o{io},filesep,o{io}],'prior','data','forward');
        prior_mul{io}=m.prior;
        data_mul{io}=m.data;
        forward_mul{io}=m.forward;
        
    end
    
end

FS=6;
prior{1}.cax=[.12 .16];

%% PLOT MULTIPE KERNELS
im=1;
id=1;
N=length(length(f_mul{im}.im));
figure(20);clf
colormap(sippi_colormap(4))

for j=1:length(f_mul{im}.im)
    j
    subplot(1,10,j)
    
    m0=m_ref{1}.*0.*NaN;
    iim=f_mul{im}.im{j};
    m0(iim)=m_ref{1}(iim);
    %imagesc(prior{1}.x,prior{1}.y,m_ref{1});
    im_handle=imagesc(prior{1}.x,prior{1}.y,m0);
    colormap_nan(im_handle)
    
    
    set(gca,'FontSize',6);
    axis image
    hold on
    %plot(f_mul{im}.xx(f_mul{im}.im{j}),f_mul{im}.yy(f_mul{im}.im{j}),'g.');
    plot([forward.sources(f_mul{im}.id{j},1),f_mul{im}.receivers(f_mul{im}.id{j},1)]',[f_mul{im}.sources(f_mul{im}.id{j},2),f_mul{im}.receivers(f_mul{im}.id{j},2)]','k-*','LineWidth',.1,'MarkerSize',.1)
    hold off
    drawnow
    set(gca,'FontSize',FS-2);
    axis image;caxis(prior{1}.cax)
    t(i+1+N)=title(sprintf('%s) subset %d',char(96)+j,j),'FontSize',FS);
end

%print_mul([txt,'_subsets'],5,0,0,600);

%% PLOT 1D MARGINAL OF MODELING ERROR
figure(15);clf;
hx=linspace(-4,4,1001);
hx=linspace(-2,2,1001);
clear hc L;
for i=1:length(dd);
    hc(i,:)=hist(dd{i}(:),hx);
    L{i}=tit{i};
end
hc=hc./prod(size(dd{1}));
hc(:,1)=NaN;
hc(:,end)=NaN;
plot(hx,hc,'LineWidth',2);

legend(L)
xlabel('Data residual (ns)')
ylabel('Frequency')
ppp(12,8,FS+2,2,2)
print_mul(sprintf('%s_1Dmarg',txt),5,0,0,600);
%% PLOT POSTERIOR REALIZATIONS
n_plot=8;

x=prior{1}.x;nx=length(x);
y=prior{1}.y;ny=length(y);
nr=length(r_all);;
figure(22);set_paper('portrait');clf;

nx_s=max([10 n_plot])
ny_s=max([nr 8])

for i=1:nr
    ns=size(r_all{i},1);
    %i_plot=ceil(linspace(1,ns,2*n_plot));i_plot=i_plot((n_plot+1):end);
    i_plot=ceil(linspace(ceil(ns/20),ns,n_plot));
    
    for j=1:n_plot;
        subplot(ny_s,nx_s,(i-1)*nx_s+j);
        imagesc(x,y,reshape(r_all{i}(i_plot(j),:),ny,nx));
        axis image;
        set(gca,'FontSize',FS-2);
        axis image;caxis(prior{1}.cax)
        %t=title(sprintf('# %d',j));
        text(max(x),0,sprintf('# %d',j),'FontSize',FS,'HorizontalAlignment','right','VerticalAlignment','bottom')
        
        if j==1;
            t=text(-6,-2.5,sprintf('%s) g_{%s}',char(96+i),tit{i}));
            set(t,'HorizontalAlignment','left','FontSize',FS+2,'FontWeight','Bold');
        end
    end
    drawnow;
end
print_mul([txt,'_PosteriorReals'],1,0,0,600);

%% FORWARD CPU TIME
for io=1:length(o);
    [d,~]=sippi_forward(m_ref,forward_mul{io},prior_mul{io});
end
for io=1:length(o);
    disp(io)
    nit=100;
    t1=now;
    for i=1:nit;
        %m=sippi_prior(prior_mul{io});
        [d,~]=sippi_forward(m_ref,forward_mul{io},prior_mul{io});
    end
    t2=now;
    T_elapsed(io)=((t2-t1)*3600*24)/nit;
    
end
nit2=5;
t1=now;
for i=1:nit2;
    disp(i)
    [d,~]=sippi_forward(m_ref,forward,prior);
end
t2=now;
T_elapsed(io+1)=((t2-t1)*3600*24)/nit2;



%% MEAN AND VARIANCE
figure(1);clf;
set(gca,'FontSize',FS-4);
N=length(tit);
cax_var=[0 0.2].*1e-3;
cax_std=[0 0.008];%sqrt(cax_var);

NR=4;
subplot(NR,N+1,1);
imagesc(prior{1}.x,prior{1}.y,m_ref{1});
set(gca,'FontSize',FS-2);
axis image;caxis(prior{1}.cax)
t(1)=title(sprintf('%s) m_{ref}',char(96)+1),'FontSize',FS);

for i=1:N
    subplot(NR,N+1,i+1);
    imagesc(prior{1}.x,prior{1}.y,em{i});
    set(gca,'FontSize',FS-2);
    axis image;caxis(prior{1}.cax)
    t(i+1)=title(sprintf('%s) g_{%s}',char(96)+i+1,tit{i}),'FontSize',FS);
    if i==N; colorbar_shift;end
    
    subplot(NR,N+1,i+1+(N+1));
    imagesc(prior{1}.x,prior{1}.y,sqrt(ev{i}));
    set(gca,'FontSize',FS-2);
    axis image;caxis(cax_std)
    t(i+1+N)=title(sprintf('%s) g_{%s}',char(96)+i+1+N,tit{i}),'FontSize',FS);
    if i==N; colorbar_shift;end
    
end
for i=1:length(t);
    set(t(i),'HorizontalAlignment','left','Position',[-1 -.4],'FontSize',FS);
end
colormap(sippi_colormap(4))

print_mul([txt,'_CompareEtypes']);

%% Figure 2
figure(2);clf;
c=[-1 1].*.05;
for  i=1:length(data_mul)
    
    mean_std(i)=mean(sqrt(diag(data_mul{i}{1}.Ct)));
    
    subplot(2,N,i);
    imagesc(data_mul{i}{1}.Ct);
    axis image
    caxis(c)
    %if i==1; c=caxis;else; caxis(c);end
    set(gca,'FontSize',FS-2);
    t=title(sprintf('%s) g_{%s}',char(96+i),tit{i}),'FontSize',FS);
    set(t,'HorizontalAlignment','left','Position',[-1 -.4],'FontSize',FS);
    if i==N; colorbar_shift;end
    
end
colormap(jet);
print_mul([txt,'_CompareCt']);

disp('Mean STD of the modeling error [ns]')
disp(sprintf('  %4.2f',mean_std));

%% Mean Error
i_skip=10;
for i=1:length(o)
    Nr=size(r_all{i},1);
    mm=repmat(m_ref{1}(:)',Nr-i_skip,1);
    md=mm-r_all{i}((1+i_skip):end,:);
    MAE(i)=mean(abs(md(:)));
    MSTD(i)=std(md(:));
    
    dt_mean(i) = mean(dt{i})
    Ct_mean(i) = mean(diag(Ct{i}))
    
end

figure(5);clf,
bar(MSTD);
xlabel('Type')
ylabel('\sigma(\Delta m) (m/ns)')
ylim([.7 1].*0.013)
set(gca,'xtick',1:1:length(tit),'xticklabel',tit);
ppp(8,4,FS,2,1)
grid on
print_mul([txt,'_Compare_MSTD']);

figure(6);clf,
bar(MAE);
xlabel('Type')
ylabel('Mean absolute error (m/ns)')
ylim([min(MAE)*.95 max(MAE)*1.05])
set(gca,'xtick',1:1:length(tit),'xticklabel',tit);
ppp(8,4,FS,2,1)
grid on
print_mul([txt,'_Compare_MAE']);

%% Mean posterior variance
figure(7);clf,
bar(MEAN_POST_VAR);
set(gca,'FontSize',FS);
xlabel('Type','FontSize',FS+2)
ylabel('Mean posterior variance')
set(gca,'xtick',1:1:length(tit),'xticklabel',tit);
ppp(8,4,FS,2,1)
grid on
print_mul([txt,'_Compare_PVA']);





%% CC
clear cc_j;
for i=1:length(o)
    Nr=size(r_all{i},1);
    i_skip=20;
    try
        for j=1:(Nr-i_skip);
            mj=r_all{i}(j+i_skip,:);
            m2=m_ref{1};
            cc=corrcoef([mj(:),m2(:)]);
            cc_j(i,j)=cc(2);
        end
        MCC(i)=mean(cc_j(i,1:j));
    end
end
MCC
figure(8);
bar(MCC);
set(gca,'FontSize',FS);
xlabel('Type','FontSize',FS+2)
ylabel('Mean posterior correlation coefficient')
set(gca,'xtick',1:1:length(tit),'xticklabel',tit);
ylim([0.4 0.8 ])
ppp(8,4,FS,2,1)
grid on
print_mul([txt,'_Compare_MCC']);


%% SHOW TABLE



txt1=sprintf('%10s&','-');
for i=1:length(o)
    txt1=[txt1,sprintf('%10s&',tit{i})];
end
txt_tit=sprintf('%10s&','CC');
txt_cc1=sprintf('%10.2g&',MCC);
txt_cc=[txt_tit,txt_cc1];

txt_tit=sprintf('%10s&','MSTD (ns)');
txt_std1=sprintf('%10.2g&',mean_std);
txt_std=[txt_tit,txt_std1];


txt_tit=sprintf('%10s&','MAE (m/ns)');
txt_mae1=sprintf('%10.2g&',MAE);
txt_mae=[txt_tit,txt_mae1];

txt_tit=sprintf('%10s&','TIME (ns)');
txt_t1=sprintf('%10.1f&',1000*T_elapsed);
txt_t=[txt_tit,txt_t1];

txt_tit=sprintf('%10s&','TIME (SU)');
txt_t1=sprintf('%10d&',round(T_elapsed(end)./T_elapsed));
txt_su=[txt_tit,txt_t1];

txt_tit=sprintf('%10s&','MEAN DT');
txt_t1=sprintf('%10.2f&',dt_mean);
txt_dt=[txt_tit,txt_t1];

txt_tit=sprintf('%10s&','sqrt(Ct)');
txt_t1=sprintf('%10.2f&',sqrt(Ct_mean));
txt_ct=[txt_tit,txt_t1];






hr='_____________________________________________________________________________________________________';
disp(hr)
disp(txt1)
disp(hr)
disp(txt_dt)
disp(txt_ct)
disp(txt_std)
disp(txt_cc)
disp(txt_mae)
disp(txt_t)
disp(txt_su)

fid=fopen(sprintf('%s_Compare.txt',txt),'w');
fprintf(fid,'%s\n',hr);
fprintf(fid,'%s\n',txt1);
fprintf(fid,'%s\n',hr);
fprintf(fid,'%s\n',txt_dt);
fprintf(fid,'%s\n',txt_ct);

fprintf(fid,'%s\n',txt_std);
fprintf(fid,'%s\n',txt_cc);
fprintf(fid,'%s\n',txt_mae);
fprintf(fid,'%s\n',txt_t);
fprintf(fid,'%s\n',txt_su);
fprintf(fid,'%s\n',hr);

fclose(fid)




%% PLOT JACOBIAN
caxJ=[-1 1].*3;
for iray=80;%80:10:700;%320
    figure(21);clf;
    cmap=cmap_linear([1 0 0; 1 1 1;0 0 0]);
    colormap(cmap);
    for ifor=1:length(f_mul);
        ifor
        subplot(1,8,ifor);
        J{ifor}=sippi_forward_jacobian(m_ref,f_mul{ifor},prior,.0001,iray);
        J{ifor}(find(J{ifor}==0))=NaN;
        im_handle=imagesc(prior{1}.x,prior{1}.y,J{ifor});
        set(gca,'FontSize',FS-2);
        colormap_nan(im_handle)
        axis image;
        
        caxis(caxJ);
        t(ifor+1+N)=title(sprintf('%s) J(%s)',char(96)+ifor,tit{ifor}),'FontSize',FS);
        drawnow;
    end
    
    ifor=ifor+1;
    subplot(1,8,ifor);
    J{ifor}=sippi_forward_jacobian(m_ref,forward,prior,.0001,iray);
    J{ifor}(find(J{ifor}==0))=NaN;
    im_handle=imagesc(prior{1}.x,prior{1}.y,J{ifor});
    colormap_nan(im_handle)
    set(gca,'FontSize',FS-2);
    axis image;
    
    caxis(caxJ);
    t(ifor+1+N)=title(sprintf('%s) J(%s)',char(96)+ifor,'g_{fd}'),'FontSize',FS);
    
    print_mul(sprintf('%s_J%03d',txt,iray),5,0,0,600);
end

