% prior_reals_mps: Sampling MPS type prior models
rng('default');rng(6);
clear all;

nr=200;

n_multiple_grids=5;
n_cond_1d=9;
n_cond=n_cond_1d^2;

TI=channels;
TI=TI(2:2:end,2:2:end);

x=[0:1:80]*.5;
y=[0:1:40]*.5;
x=[0:1:80];
y=[0:1:40];

x=[0:2:80];
y=[0:2:40];

z=0;

n_hard=20;
n_soft=20;

use_hard=1;
use_soft=1;


%% CONDITIONAL DATA

% HARD DATA
clear d_c
ix=ceil(rand(1,n_hard).*length(x));
iy=ceil(rand(1,n_hard).*length(y));
iz=iy.*0+1;
for i=1:n_hard
    d_c(i)=TI(iy(i),ix(i));
    %d_c(i)=1;imagesc(m{1})
end
f_cond=sprintf('f_cond_%d.dat',n_hard);
d_cond=[x(ix)' y(iy)' z(iz)' d_c(:)];
write_eas(f_cond,d_cond);


% SOFT
clear d_s
ix=ceil(rand(1,n_soft).*length(x));
iy=ceil(rand(1,n_soft).*length(y));
iz=iy.*0+1;
%d_x=ceil(rand(1,n_soft).*(max(x)-min(x))+min(x));
%d_y=ceil(rand(1,n_soft).*(max(y)-min(y))+min(y));
%d_z=d_y.*0+z(1);
for i=1:n_soft
    p=rand(1)*.01;
    v=TI(iy(i),ix(i));
    if v==0;
        d_s(i)=1-p;
    else
        d_s(i)=p;
    end
end
f_soft=sprintf('f_soft_%d.dat',n_soft);
d_soft=[x(ix)' y(iy)' z(iz)' d_s(:) 1-d_s(:)];
write_eas(f_soft,d_soft);



%% SETUP DIFFERENT TYPE OF PRIOR



ip=1; 
prior{ip}.type='mps'; % MPS type
prior{ip}.method='mps_snesim_tree'; 
prior{ip}.x=x; % X array 
prior{ip}.y=y; % Y array 
prior{ip}.ti=TI;
if (use_hard), prior{ip}.hard_data=d_cond; end
if (use_soft), prior{ip}.soft_data=d_soft; end
prior{ip}.MPS.template_size=[n_cond_1d n_cond_1d 1];
prior{ip}.MPS.n_multiple_grids=n_multiple_grids;
prior{ip}.MPS.shuffle_simulation_grid=2;

ip=ip+1; 
prior{ip}.type='mps'; % MPS type
prior{ip}.method='mps_enesim'; 
prior{ip}.x=x; % X array 
prior{ip}.y=y; % Y array 
prior{ip}.ti=TI;
if (use_hard), prior{ip}.hard_data=d_cond; end
if (use_soft), prior{ip}.soft_data=d_soft; end
prior{ip}.MPS.n_max_cpdf_count=10;
prior{ip}.MPS.n_max_ite=100000;
prior{ip}.MPS.n_cond=n_cond;
prior{ip}.MPS.shuffle_simulation_grid=2;

%ip=ip+1; 
%prior{ip}.type='snesim_std'; % MPS type
%prior{ip}.x=x; % X array 
%prior{ip}.y=y; % Y array 
%prior{ip}.ti=TI;

%ip=ip+1; 
%prior{ip}.type='snesim'; % MPS type
%prior{ip}.x=x; % X array 
%prior{ip}.y=y; % Y array 
%prior{ip}.ti=TI;
%prior{ip}.hard_data=d_cond;
%[m,prior]=sippi_prior(prior);
%prior{ip}.S.nmulgrids=n_multiple_grids;
%prior{ip}.S.max_cond=n_cond;











%%
figure(1);clf;set_paper;
np=length(prior);

nr_show=3;
for ir=1:nr;
[m,prior]=sippi_prior(prior);
    
for ip=1:np;
    
    if ir<=nr_show;
        subplot(nr_show+2,np,(ir-1)*np+ip)
        imagesc(x,y,m{ip});axis image
        colormap(cmap_linear([1 0 0;1 1 1;0 0 0]))
        ylabel(sprintf('t=%4.2gs',prior{ip}.time))
        
        
        hold on
        plot(d_cond(:,1),d_cond(:,2),'w.','MarkerSize',25)
        scatter(d_cond(:,1),d_cond(:,2),15,d_cond(:,4),'filled')
        plot(d_soft(:,1),d_soft(:,2),'y.','MarkerSize',25)
        scatter(d_soft(:,1),d_soft(:,2),10,1-d_soft(:,4),'filled')
        hold off
        
        if ir==1;99
            try
                title(prior{ip}.type,'interpreter','none')
            catch
                title(prior{ip}.method,'interpreter','none')
            end
        end
        drawnow;
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
    plot(d_cond(:,1),d_cond(:,2),'w.','MarkerSize',25)
    scatter(d_cond(:,1),d_cond(:,2),15,d_cond(:,4),'filled')
    plot(d_soft(:,1),d_soft(:,2),'y.','MarkerSize',25)
    scatter(d_soft(:,1),d_soft(:,2),10,1-d_soft(:,4),'filled')
    hold off
    title(sprintf('#=%03d',ir))    
    drawnow
    
    
end
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
    title(sprintf('t=%3.2g',prior{1}.time))
end
colormap(sippi_colormap(1));
colorbar_shift;
print_mul(sprintf('prior_reals_%ss',prior{1}.type));
s=suptitle(sprintf('Independent realizations from a %s (%s) type prior',upper(prior{1}.type),prior{1}.method))
set(s,'interpreter','none')

%% CONDITIONAL TO HARD DATA

%% CONDITIONAL TO SOFT DATA

%% CONDITIONAL TO HARD AND SOFT DATA

