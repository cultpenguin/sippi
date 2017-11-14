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
% NOTE: Only works to high dimensional prior models....
%
% See also sippi_metropolis, sippi_metropolis_mulruns
%
function [cc_mix]=sippi_plot_posterior_mixing(o,txt,Nplot);

%% CHECK INPUT
Nc=length(o);
if nargin<2,
    txt=sprintf('test_mixing');
    txt=[o{1}.txt,'_mixing']
end
if nargin<3,
    Nplot=Nc;
end
io_arr=1:Nplot;

%% LOAD DATA
skip_seq_gibbs=0;
disp(sprintf('%s: reading data ...',mfilename))
for io=1:Nc;
    [reals{io},etype_mean{io},etype_var{io},reals_all{io},reals_ite{io}]=sippi_get_sample(o{io}.txt,1,15,skip_seq_gibbs);
end
load(sprintf('%s%s%s',o{1}.txt,filesep,o{1}.txt),'options');
options.mcmc.null='';
load(sprintf('%s%s%s',o{1}.txt,filesep,o{1}.txt),'prior');

if isfield(options.mcmc,'m_ref')
    is_m_ref=1;
else
    is_m_ref=0;
end
    


%% IMAGE ETYPE + REF
figure(1);clf;
j=0;
Ny=ceil((Nc+is_m_ref)/5);
Nx=ceil((Nc+is_m_ref)/Ny);
if (is_m_ref)    
    N=length(o)+1;
    j=j+1;subplot(Ny,Nx,1);
    imagesc(prior{1}.x,prior{1}.y,options.mcmc.m_ref{1});
    axis image
    caxis(prior{1}.cax)
    title('Ref')
end

for io=1:Nc;
    j=j+1;subplot(Ny,Nx,j);
    imagesc(prior{1}.x,prior{1}.y,etype_mean{io});
    axis image
    caxis(prior{1}.cax)    
    title(sprintf('Chain #%02d',io))
end
axis image
print_mul(sprintf('%s_etype',txt))   
grid on

%% logL
figure(3);clf;
for io=1:Nc;
    load([o{io}.txt,filesep,o{io}.txt,'.mat'],'C')    
    h(io)=plot(C{1}.mcmc.logL,'-');
    %sippi_plot_loglikelihood(C{1}.mcmc.logL)
    hold on
    L{io}=sprintf('Chain %d',io);
end
hold off
legend(L,'Location','NorthEastOutside')
xlabel('Iteration number')
ylabel('log-likelihood')
print_mul(sprintf('%s_logL',txt))   


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
        if (is_m_ref);
            cc=corrcoef(options.mcmc.m_ref{1}(:),reals_all{io}(i,:)');
            cc_ref{io}(j)=cc(2);
        end
        
        % i_use=ceil(nr/3);
        i_use=nr;
        for io_end=1:Nc
            cc=corrcoef(reals_all{io_end}(i_use,:),reals_all{io}(i,:)');
            cc_mix{io}{io_end}(j)=cc(2);
        end
    end
    
    for io_end=1:Nc;
        if io_end==io; 
            lw=3;
            col='r';
        else
            lw=1;
            col='k';
        end
        figure_focus(10+io);subplot(2,1,1);
    
        p2(io_end)=plot(i_ax,cc_mix{io}{io_end},'k-','LineWidth',lw);
        set(p2(io_end),'color',col);
        ih=ih+1;leg{ih}=sprintf('CC_{%d,%d}',io,io_end);
        if ih==1; hold on; end
        
        hx=linspace(-1,1,41);
        [h(io_end,:)]=hist(cc_mix{io}{io_end},hx);
        
    end
    p_all=p2;
    if (is_m_ref)
        figure_focus(10+io);subplot(2,1,1);
    
        p=plot(i_ax,cc_ref{io},'r-','LineWidth',2);
        ih=ih+1;leg{ih}=sprintf('CC_{%d,ref}',io);
        if ih==1; hold on; end
        
        [h(Nc+1,:)]=hist(cc_ref{io},hx);
        
        p_all= [p_all p];
        
        
    end


    figure_focus(10+io);
    subplot(2,1,1);
    hold off
    legend([p_all],leg,'Location','NorthEastOutside')
    title(sprintf('%s - Chain %02d vs others',txt,io),'Interpreter','None')
    ylabel(sprintf('Correlation coefficient to Last model of chain %02d',io))
    xlabel('Iteration number')
    
    arr1=io;
    arr2=[setxor(1:Nc,io)];
    %if (is_m_ref)
    %    arr2=[arr2,Nc+1]
    %end
    
    subplot(2,1,2);
    %    figure_focus(200+io);
    b1=bar(hx,h(arr1,:),'k');
    hold on
    b2=plot(hx,h(arr2,:)');
    if (is_m_ref)
        b3=bar(hx,h(Nc,:)','w');   
        set(get(b3,'Children'),'FaceAlpha',0.3)
        legend([b1;b2;b3],leg{[arr1,arr2,Nc+1]},'Location','NorthEastOutside')
    else
        legend([b1;b2],leg{[arr1,arr2]},'Location','NorthEastOutside')
    end
    hold off
    xlabel(sprintf('Correlation coefficient\n(Last model of chain %02d to posterior realization of all other chains)',io))
    ylabel('Relative frequency')
    drawnow;
    print_mul(sprintf('%s_mixing_C%02d',txt,io))
    
    
end


