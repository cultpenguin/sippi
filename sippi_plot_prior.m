function sippi_plot_prior(prior,im_arr,n_reals,caxis,supt);
% sippi_plot_prior : generate a sample of the prior
%
% Call :
%    sippi_plot_prior(prior,ip,n_reals,cax,supt);
%
%
%  See also: sippi_plot_posterior, sippi_plot_model
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
        if prior{im_arr(j)}.ndim<2
            n_reals(j)=10000;
        else
            n_reals(j)=15;
        end
    end
end

if ~exist('cax','var');
    try
        cax=prior{1}.cax;
    end
end


if length(n_reals)==1;
    n_reals=ones(1,length(im_arr)).*n_reals;
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


for im=im_arr;
    % CLEAR FIGURES
    figure_focus(99+im);clf;
    set_paper('landscape');
    
    if prior{im}.ndim<2
        % HISTOGRAM
        p{1}=prior{im};
        m_reals=zeros(1,n_reals(im));
        for i=1:n_reals(im);            
            m=sippi_prior(p);
            m_reals(i)=m{1};       
        end

        hist(m_reals,30);
        xlabel(prior{im}.name)
        
        if isfield(prior{im},'mlim');
            set(gca,'xlim',prior{im}.mlim);
        end
        title(sprintf('m#%d',im))
        
    else
        %% SUBPLOTS       
        
        for i=1:n_reals(im);
            progress_txt(i,n_reals,'generating prior sample')
            if (prior{im}.dim(1)>max(prior{im}.dim(2:3)))
                nsp_y=5;
                nsp_x=ceil(n_reals(im)/nsp_y);
            else
                nsp_x=5;
                nsp_y=ceil(n_reals(im)/nsp_x);
            end
            
            try;prior{1}.seed=prior{1}.seed+1;end;
            m_prior=sippi_prior(prior);
            figure_focus(99+im);set_paper('portrait');%set_paper('landscape');
            subplot(nsp_y,nsp_x,i);
            use_colorbar=0;
            
            if ((n_reals==i)&(i==(1*nsp_x)))|(i==(2*nsp_x));
                use_colorbar=1;
            end
            
            sippi_plot_model(prior,m_prior,1:nm,use_colorbar);
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