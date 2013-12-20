function sippi_plot_prior(prior,im_arr,n_reals,caxis,supt);
% sippi_plot_prior Plot a sample of the prior in SIPPI
%
% Call :
%    sippi_plot_prior(prior,ip,n_reals,cax,supt);
%
%  See also sippi_plot_posterior, sippi_plot_model
%

cwd=pwd;

%% DATA
if isstr(prior)
    % load prior from MATFILE
    fname=prior;
    try
        load([fname,'.mat']);
    catch
        cd(fname);
        load([fname,'.mat']);
    end
end
prior=sippi_prior_init(prior);

if nargin>1
    if isempty(im_arr)
        im_arr=1:length(prior);
    end
else    
    im_arr=1:1:length(prior);
end
   
if nargin<5,
    supt=0;
end

if ~exist('n_reals','var');
    for j=1:length(im_arr);
        if prior{im_arr(j)}.ndim<1
            %n_reals(j)=1000;
            n_reals(im_arr(j))=1000;
        elseif prior{im_arr(j)}.ndim<2
            %n_reals(j)=30;
            n_reals(im_arr(j))=30;
        else
            %n_reals(j)=15;
            n_reals(im_arr(j))=15;
        end
    end
end

if length(n_reals)==1;
    %n_reals=ones(1,length(im_arr)).*n_reals;
    n_reals=ones(1,length(prior)).*n_reals;
end
   
if ~exist('cax','var');
    try
        cax=prior{1}.cax;
    end
end

 

nm=length(prior);
j=0;
for im=im_arr
    j=j+1;
    try
        prior{im}.fftma_options=rmfield(prior{im}.fftma_options,'z_rand');
    end
    prior{im}.ndim=sum(find(prior{im}.dim>1));
end

prior_master=zeros(1,length(prior));
for im=1:length(prior);
    if isfield(prior{im},'prior_master');
        prior_master(prior{im}.prior_master)=1;
    end
end

for im=im_arr;
    % CLEAR FIGURES
    f_id=90+im;
    figure_focus(f_id);clf;
    set_paper('landscape');
    
    
    % FIND SCALE/ORIENTATION
    ax_lscape=1;
    try
        % if 'lim' is set
        if prior{im}.lim(1)<max(prior{im}.lim(2:3))
            ax_lscape=0;
        end
    end
    try
        % if 'daspect' is set
        r=prior{im}.lim./prior{im}.daspect;
        if r(1)<max(r(2:3))
            ax_lscape=0;
        end
    end
    
    if isfield(prior{im},'name');
        title_txt=sprintf('m%d: %s',im,prior{im}.name);
    else
        title_txt=sprintf('m%d',im);
    end
    
    
    if prior{im}.ndim<1
        % SCALAR --> HISTOGRAM
        clear p;
        p{1}=prior{im};
        m_reals=zeros(1,n_reals(im));
        for i=1:n_reals(im);
            m=sippi_prior(p);
            m_reals(i)=m{1};       
        end
        
        hist(m_reals,30);
        ylim=get(gca,'ylim');
        hold on;
        try
            % requires statistics toolbox
            plot([1 1].*quantile(m_reals,.025),ylim,'r--','linewidth',2);
            plot([1 1].*quantile(m_reals,.5),ylim,'r-','linewidth',2);
            plot([1 1].*quantile(m_reals,.975),ylim,'r--','linewidth',2);
        end
        hold off
        xlabel(prior{im}.name)
        if isfield(prior{im},'cax');
            set(gca,'xlim',prior{im}.cax);
        end
               
        title(title_txt)
        
    elseif prior{im}.ndim<2
        % 1D
        clear p;clear m_reals;
        if (prior_master(im)==1);
            % FULL PRIOR REALS
            for i=1:n_reals(im);
                m=sippi_prior(prior);
                m_reals(i,:)=m{im};
            end
        
        else
            % ONLY REALS FROM ONE/CURRENT PRIOR
            p{1}=prior{im};
            for i=1:n_reals(im);
                m=sippi_prior(p);
                m_reals(i,:)=m{1};
            end
        end
        
        plot(prior{im}.x,m_reals,'k-');
        hold on
        plot(prior{im}.x,quantile(m_reals,.025),'r--','linewidth',2);
        plot(prior{im}.x,quantile(m_reals,.5),'r-','linewidth',2);
        plot(prior{im}.x,quantile(m_reals,.975),'r--','linewidth',2);
        hold off
        
        if isfield(prior{im},'cax');
            set(gca,'ylim',prior{im}.cax);
        end
        
        xlabel('X')
        ylabel(prior{im}.name)
        title(title_txt)
    else
        %% SUBPLOTS       
        for i=1:n_reals(im);
            progress_txt(i,n_reals(im),'generating prior sample')
             if ax_lscape==1;
                nsp_y=5;
                nsp_x=ceil(n_reals(im)/nsp_y);
            else
                nsp_x=5;
                nsp_y=ceil(n_reals(im)/nsp_x);
             end
            
            clear p;clear m_reals;
            try;clear m_prior;end  
            if (prior_master(im)==1);
                % FULL PRIOR REALS
                m_prior_full=sippi_prior(prior);
                p{1}=prior{im};
                m_prior{1}=m_prior_full{im};
            else
                p{1}=prior{im};
                try;p{1}.seed=p{1}.seed+1;end;
                m_prior=sippi_prior(p);
            end
            figure_focus(f_id);
            
            subplot(nsp_y,nsp_x,i);
            
            use_colorbar=0;
            i_cb=ceil((nsp_y+1)/2)*nsp_x;
            if i==i_cb; use_colorbar=1;end
            
            sippi_plot_model(p,m_prior,1,use_colorbar,f_id);
        end
        %%%
        
    end
    
    if supt==1,
        sp=suptitle(sprintf('m%d : %s',im,prior{im}.name));
        set(sp,'interp','none')
    end
    try
        print_mul(sprintf('%s_m%d_prior_sample',fname,im))
    catch
        print_mul(sprintf('m%d_prior_sample',im))
    end

    
end

cd(cwd);

return



%    end


x=forward.x;
y=forward.y;
z=forward.z;


%% PLOT NOISE MODEL
%try;
figure(21);set_paper('landscape');clf;

i_use=forward.i_use;

ha=tight_subplot(1,2,[.08 .08],[.1 .1],[.1 .1]);

%subplot(1,2,1);
axes(ha(2));
ax1=gca;
set(ax1,'FontSize',16);
imagesc(forward.CD(i_use,i_use));
xlabel('Data #')
ylabel('Data #')
%axis image
cb=colorbar;
pos=get(ax1,'position');
set(get(cb,'Ylabel'),'String','Covariance')
title('CD')
axis image

%subplot(1,2,2);
axes(ha(1));
ax2=gca;
set(ax2,'FontSize',16);
plot(forward.dt0(i_use));
xlabel('Data #')
ylabel('d_0 (ns)')
set(ax2,'xlim',[0 length(i_use)])
pos2=get(ax2,'position');
pos = plotboxpos(ax1);
pos2(2)=pos(2);
pos2(4)=pos(4);
set(ax2,'position',pos2);
%end
print_mul(sprintf('%s_noisemodel',fname))

return
%% PLOT REALS
figure(12);set_paper('landscape');clf;

nsp_x=5;
nsp_y=ceil(n_reals/nsp_x);
for i=1:n_reals
    subplot(nsp_y,nsp_x,i);
    imagesc(x,y,m_prior(:,:,i));
    caxis(cax)
    if (i==(2*nsp_x));
        colorbar_shift;
    end
    axis image;

end
colormap(1-gray)
try
    print_mul(sprintf('%s_prior_sample',fname))
end