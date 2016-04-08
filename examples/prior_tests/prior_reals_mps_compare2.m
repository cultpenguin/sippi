% prior_reals_mps: Sampling MPS type prior models
rng('default');rng(2);
clear all;

%% BASIC SETTINGS
TI=channels;
TI=TI(2:2:end,2:2:end);

x=[0:1:80]*.5;y=[0:1:40]*.5;
x=[0:2:80];y=[0:2:40];
x=[0:1:80];y=[0:1:40];
%x=[0:3:240];y=[0:3:120];

nr=100;
%nr=3;

use_hard=1;n_hard=10;
use_soft=0;
use_soft_random=1;n_soft=50;
use_soft_random=0;n_soft=length(x)*length(y);

shuffle_simulation_grid=2;
n_multiple_grids=3;
n_cond_1d=7;
n_cond=n_cond_1d^2;
z=0;

cmap=(cmap_linear([1 0 0;1 1 1;0 0 0]));
cax=[0 1];
%% SETUP DIFFERENT TYPE OF PRIOR

ip=1;
prior{ip}.type='mps'; % MPS type
prior{ip}.method='mps_snesim_tree';
prior{ip}.name='MPS - SNESIM_TREE'; % MPS type
prior{ip}.x=x; % X array
prior{ip}.y=y; % Y array
prior{ip}.ti=TI;
prior{ip}.MPS.template_size=[n_cond_1d n_cond_1d 1];
prior{ip}.MPS.n_multiple_grids=n_multiple_grids;
prior{ip}.MPS.shuffle_simulation_grid=shuffle_simulation_grid;

ip=ip+1;
prior{ip}.type='mps'; % MPS type
prior{ip}.method='mps_enesim_general';
prior{ip}.name='MPS - ENESIM/DSIM'; % MPS type
prior{ip}.x=x; % X array
prior{ip}.y=y; % Y array
prior{ip}.ti=TI;
%if (use_hard), prior{ip}.hard_data=d_hard; end
%if (use_soft), prior{ip}.soft_data=d_soft; end
prior{ip}.MPS.n_max_cpdf_count=1; % D SAMP
prior{ip}.MPS.n_max_ite=100000;
prior{ip}.MPS.n_cond=n_cond;
prior{ip}.MPS.shuffle_simulation_grid=shuffle_simulation_grid;

ip=ip+1;
prior{ip}.type='mps'; % MPS type
prior{ip}.method='mps_enesim';
prior{ip}.name='MPS - ENESIM'; % MPS type
prior{ip}.x=x; % X array
prior{ip}.y=y; % Y array
prior{ip}.ti=TI;
prior{ip}.MPS.n_max_cpdf_count=100; % D SAMP
prior{ip}.MPS.n_max_ite=100000;
prior{ip}.MPS.n_cond=n_cond;
prior{ip}.MPS.shuffle_simulation_grid=shuffle_simulation_grid;

ip=ip+1;
prior{ip}.type='snesim_std'; % MPS type
prior{ip}.name='SGeMS snesim'; % MPS type
prior{ip}.x=x; % X array
prior{ip}.y=y; % Y array
prior{ip}.ti=TI;
prior{ip}.S.XML.parameters.Nb_Multigrids_ADVANCED.value=n_multiple_grids+1;
prior{ip}.S.XML.parameters.Max_Cond.value=n_cond;

ip=ip+1;
prior{ip}.type='snesim'; % MPS type
prior{ip}.name='Fortran snesim'; % MPS type
prior{ip}.x=x; % X array
prior{ip}.y=y; % Y array
prior{ip}.ti=TI;
prior{ip}.S.XML.parameters.nmulgrids=n_multiple_grids+1;
prior{ip}.S.XML.parameters.max_Cond=n_cond;

%%
%[m,prior]=sippi_prior(prior);

%return
%%
i_ref=length(prior);
p_ref{1}=prior{i_ref};
p_ref{1}.MPS.rseed=3;
p_ref{1}.MPS.template_size=[9 9 1];
p_ref{1}.MPS.n_cond=81;
m_ref=sippi_prior(p_ref);

%template_size=prior{1}.MPS.template_size;
%prior{1}.MPS.rseed=1;
%prior{1}.MPS.template_size=[9 9 1];
%[m_ref,prior]=sippi_prior(prior);
%prior{1}.MPS.template_size=template_size;


%% CONDITIONAL DATA

% HARD DATA
clear d_h
ix=ceil(rand(1,n_hard).*length(x));
iy=ceil(rand(1,n_hard).*length(y));
iz=iy.*0+1;
for i=1:n_hard
    d_h(i)=m_ref{1}(iy(i),ix(i));
end
f_cond=sprintf('f_cond_%d.dat',n_hard);
d_hard=[x(ix)' y(iy)' z(iz)' d_h(:)];
write_eas(f_cond,d_hard);


% SOFT
clear d_s
if use_soft_random==1
    ix=ceil(rand(1,n_soft).*length(x));
    iy=ceil(rand(1,n_soft).*length(y));
    iz=iy.*0+1;
else
    [xx,yy]=meshgrid(1:length(x),1:length(y));
    ix=xx(1:n_soft);
    iy=yy(1:n_soft);
    iz=ix.*0+1;
    %d_x=ceil(rand(1,n_soft).*(max(x)-min(x))+min(x));
    %d_y=ceil(rand(1,n_soft).*(max(y)-min(y))+min(y));
    %d_z=d_y.*0+z(1);
end
for i=1:n_soft
    
    % horizontal gradient
    fac=(1-ix(i)/length(x));;
    p=fac.*1*.5;
    
    % random
    fac=1;
    p=fac.*rand(1)*.5;
    
    % constant weight
    %fac=.70;
    %p=fac.*1*.5;
    
    
    v=m_ref{1}(iy(i),ix(i));
    if v==0;
        d_s(i)=0.5+p;
    else
        d_s(i)=0.5-p;
        %d_s(i)=p.*.1;
    end
end
f_soft=sprintf('f_soft_%d.dat',n_soft);
d_soft=[x(ix)' y(iy)' z(iz)' d_s(:) 1-d_s(:)];
write_eas(f_soft,d_soft);


%% SELECT HARD/SOFT DATA
for ip=1:length(prior);
    if (use_hard), prior{ip}.hard_data=d_hard; end
    if (use_soft), prior{ip}.soft_data=d_soft; end
end

%%
figure(2);clf;set_paper;
np=length(prior);
nr_show=3;
for ir=1:nr;
    [m,prior]=sippi_prior(prior); % WHY DOES THIS NOT WORK!!!
    %[m]=sippi_prior(prior);
    
    subplot(nr_show+2,np,(nr_show+1)*np+1);
    imagesc(x,y,m_ref{1});axis image;
    colormap(cmap);caxis(cax)
    title('Reference model')
    ax=axis;
    subplot(nr_show+2,np,(nr_show+1)*np+2);
    scatter(d_hard(:,1),d_hard(:,2),15,d_hard(:,4),'filled');axis image;
    box on
    colormap(cmap);caxis(cax)
    set(gca,'ydir','reverse')
    title('Hard data')
    axis(ax)
    subplot(nr_show+2,np,(nr_show+1)*np+3);
    scatter(d_soft(:,1),d_soft(:,2),15,d_soft(:,5),'filled');axis image;
    box on
    colormap(cmap);caxis(cax)
    set(gca,'ydir','reverse')
    title('Soft data, P(m=1)')
    axis(ax)
    colorbar_shift;
    
    for ip=1:np;
        
        if ir<=nr_show;
            figure_focus(2);
            subplot(nr_show+2,np,(ir-1)*np+ip)
            imagesc(x,y,m{ip});axis image
            colormap(cmap)
            try;ylabel(sprintf('t=%4.2gs',prior{ip}.time));end
            
            hold on
            if (use_hard)
                plot(d_hard(:,1),d_hard(:,2),'w.','MarkerSize',25)
                scatter(d_hard(:,1),d_hard(:,2),15,d_hard(:,4),'filled')
            end
            if (n_soft<30)&&(use_soft)
                plot(d_soft(:,1),d_soft(:,2),'y.','MarkerSize',25)
                scatter(d_soft(:,1),d_soft(:,2),10,1-d_soft(:,4),'filled')
            end
            hold off
            
            if ir==1;
                try
                    title(prior{ip}.name,'interpreter','none');
                catch
                    title(prior{ip}.type,'interpreter','none');
                end
            end
        end
        
        %% MEAN
        % update mean
        if ir==1;
            em{ip}=m{ip};
        else
            em{ip}=em{ip}+m{ip};
        end
        
        
        subplot(nr_show+2,np,nr_show*np+ip)
        imagesc(x,y,em{ip}./ir);axis image;caxis([0 1])
        hold on
        if (use_hard)
            plot(d_hard(:,1),d_hard(:,2),'w.','MarkerSize',25)
            scatter(d_hard(:,1),d_hard(:,2),15,d_hard(:,4),'filled')
        end
        if (n_soft<30)&&(use_soft)
            plot(d_soft(:,1),d_soft(:,2),'y.','MarkerSize',25)
            scatter(d_soft(:,1),d_soft(:,2),10,1-d_soft(:,4),'filled')
        end
        hold off
        title(sprintf('#=%03d',ir))
        drawnow
    end
    
    drawnow;
end




print_mul(sprintf('%s_dx%g_H%d_S%d',mfilename,100*(x(2)-x(1)),use_hard,use_soft))
return

%% UNCONDITIONAL
figure(10);clf
for i=1:5;
    [m,prior]=sippi_prior(prior);
    subplot(1,5,i);
    imagesc(prior{1}.x,prior{1}.y,m{1});
    axis image
    %xlabel('X')
    %ylabel('Y')
    %    title(sprintf('t=%3.2g',prior{1}.time))
end
colormap(sippi_colormap(1));
colorbar_shift;
print_mul(sprintf('prior_reals_%ss',prior{1}.type));
s=suptitle(sprintf('Independent realizations from a %s (%s) type prior',upper(prior{1}.type),prior{1}.method))
set(s,'interpreter','none')

%% CONDITIONAL TO HARD DATA

%% CONDITIONAL TO SOFT DATA

%% CONDITIONAL TO HARD AND SOFT DATA

