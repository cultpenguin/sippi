% sippi_plot_posterior_mixing: Plots to use for analysis if mxixing of seperate metropolis
% runs has takem place
%
% % Example:
% % 1. Run several different runs of sippi_metropolis using 
% [o{1}]=sippi_metropolis(data,prior,forward);
% [o{2}]=sippi_metropolis(data,prior,forward);
% [o{3}]=sippi_metropolis(data,prior,forward);
% % automatically
% options.nruns=3;
% [o]=sippi_metropolis(data,prior,forward);
% 2. call sippi_plot_posterio_mixing
% sippi_plot_posterior_mixing(o),
%
% See also sippi_metropolis, sippi_metropolis_mulruns
%
function sippi_plot_posterior_mixing(o,txt,io_arr);

%% CHECK INPUT
if nargin<2,
    txt=sprintf('test_mixing');
end
if nargin<3,
    io_arr=1:length(o);
end

%% LOAD DATA
skip_seq_gibbs=0;
for io=1:length(o);
    [reals{io},etype_mean{io},etype_var{io},reals_all{io},reals_ite{io}]=sippi_get_sample(o{io}.txt,1,15,skip_seq_gibbs);
end
%load(sprintf('%s%s%s',o{1}.txt,filesep,o{1}.txt),'options');
options.mcmc.null='';
load(sprintf('%s%s%s',o{1}.txt,filesep,o{1}.txt),'prior');

%% IMAGE ETYPE + REF
figure(2);clf;
N=length(o);
j=0;
if isfield(options.mcmc,'m_ref')
    N=length(o)+1;
    j=j+1;subplot(1,N,1);
    imagesc(prior{1}.x,prior{1}.y,options.mcmc.m_ref{1});
    axis image
    caxis(prior{1}.cax)
    title('Ref')
end

for io=1:length(o);
    j=j+1;subplot(1,N,j);
    imagesc(prior{1}.x,prior{1}.y,etype_mean{io});
    axis image
    caxis(prior{1}.cax)    
    title(sprintf('RUN #%02d',io))
end
axis image
print_mul(sprintf('%s_etype',txt))   

%% COMPUTE CC
clear c* h*
nr=size(reals_all{1},1);
i_sample=o{1}.mcmc.i_sample;
nite=o{1}.mcmc.nite;
i_ax=i_sample:i_sample:nite;
clear h* cc* leg*
for io=io_arr%;:length(o);;
    disp(sprintf('%s: Analyzing chain %02d',mfilename,io));
    clear h leg
    ih=0;
    figure(10+io);clf;set_paper('landscape')
    ii=[ceil(nr/10):1:nr];
    i_ax=i_sample:i_sample:nite;
    i_ax=i_ax(ii);
    j=0;
    for i=ii
        j=j+1;
        if isfield(options.mcmc,'m_ref')
            cc=corrcoef(options.mcmc.m_ref{1}(:),reals_all{io}(j,:)');
            cc_ref{io}(j)=cc(2);
        end
        
        % i_use=ceil(nr/3);
        i_use=nr;
        for io_end=1:length(o);;
            cc=corrcoef(reals_all{io_end}(i_use,:),reals_all{io}(j,:)');
            cc_mix{io}{io_end}(j)=cc(2);
        end
    end
    
    if isfield(options.mcmc,'m_ref')
        figure_focus(10+io);subplot(2,1,1);
    
        p=plot(i_ax,cc_ref{io},'k-','LineWidth',2);
        try;set(p,'color',col{io});end
        ih=ih+1;leg{ih}=sprintf('Ref_%d',io);
        if ih==1; hold on; end
    end
    for io_end=1:length(o);;
        if io_end==io; 
            lw=3;
            col='r';
        else
            lw=1;
            col='k';
        end
        figure_focus(10+io);subplot(2,1,1);
    
        p2=plot(i_ax,cc_mix{io}{io_end},'k-','LineWidth',lw);
        try;set(p2,'color',col);end
        ih=ih+1;leg{ih}=sprintf('mix_{%d,%d}',io,io_end);
        if ih==1; hold on; end
        
        hx=linspace(-1,1,21);
        [h(io_end,:)]=hist(cc_mix{io}{io_end},hx);
        
        
    end

    
    figure_focus(10+io);
    subplot(2,1,1);
    hold off
    legend(leg,'Location','NorthEastOutside')
    title(sprintf('%s - Chain %02d vs others',txt,io),'Interpreter','None')
    ylabel(sprintf('Correlation coefficient to Last model of chain %02d',io))
    xlabel('Iteration number')
    
    arr1=io;
    arr2=setxor(1:length(o),io);
    
    subplot(2,1,2);
    %    figure_focus(200+io);
    b1=bar(hx,h(arr1,:),'k');
    hold on
    b2=plot(hx,h(arr2,:)');
    hold off
    legend([b1;b2],leg{[arr1,arr2]},'Location','NorthEastOutside')
    xlabel(sprintf('Correlation coefficient\n(Last model of chain %02d to posterior realization of all other chains)',io))
    ylabel('Relative frequency')
    drawnow;
    print_mul(sprintf('%s_mixing_C%02d',txt,io))
    
    
end


