clear all;
io=0;
%io=io+1;O{io}.txt='sm_20181204_2135_AM13_P1';
%io=io+1;O{io}.txt='sm_20181204_2322_AM13_P2';
%io=io+1;O{io}.txt='sm_20181205_0110_AM13_P3';
%io=io+1;O{io}.txt='sm_20181205_0412_AM13_P4';
io=io+1;O{io}.txt='sm_20181205_0614_AM13_P5';
%io=io+1;O{io}.txt='sm_20181205_1046_AM13_P6';
cax=[.12 .18];
FS=6;
for io=1:length(O);
    mat=[O{io}.txt,filesep,O{io}.txt];
    load(mat,'prior')    
    load(mat,'forward')    
    load(mat,'data')    
    
    Nshow=5;
    N=7;
    figure_focus(1);clf;
    for i=1:Nshow;
        m=sippi_prior(prior);
        figure_focus(1);
        subplot(3,N,N+i);
        imagesc(prior{1}.x,prior{1}.y,m{1})
        axis image;caxis(cax);set(gca,'FontSize',FS);
    end
    colorbar_shift;
    print_mul(sprintf('%s_prior',O{io}.txt))
    
    [reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(O{io}.txt,1,5,1);
    figure_focus(2);clf;
    for i=1:Nshow;
        m=sippi_prior(prior);
        figure_focus(2);
        subplot(3,N,N+i);
        imagesc(prior{1}.x,prior{1}.y,reals(:,:,i))
        axis image;caxis(cax);set(gca,'FontSize',FS);
    end
    colorbar_shift;
    print_mul(sprintf('%s_post',O{io}.txt))
    
    %%
    figure_focus(3);clf;
    subplot(3,N,N+1);
    imagesc(prior{1}.x,prior{1}.y,etype_mean)
    axis image;caxis(cax);set(gca,'FontSize',FS);
    subplot(3,N,N+2);
    imagesc(prior{1}.x,prior{1}.y,sqrt(etype_var))    
    axis image;caxis([0 0.02]);set(gca,'FontSize',FS);
    colorbar_shift;
    print_mul(sprintf('%s_stat',O{io}.txt))
    
        
    
    
end