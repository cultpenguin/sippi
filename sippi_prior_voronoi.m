% sippi_prior_voronoi: 
%
% TMH/2014
%
% Ex: 
%   prior{1}.type='voronoi';
%   prior{1}.x=1:1:20;
%   prior{1}.y=1:1:20;
%   prior{1}.cells_N_min=2;
%   prior{1}.cells_N_max=100;
%   prior{1}.cells_N=10;
%   [m,prior]=sippi_prior(prior);
%   sippi_plot_prior(prior,m);
%
% See also: sippi_prior_init, sippi_prior
%
function [m_propose,prior]=sippi_prior_voronoi(prior,m_current,ip);
if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

% number of voronoi cells
if ~isfield(prior{ip},'cells_N_min');prior{ip}.cells_N_min=10;end
if ~isfield(prior{ip},'cells_N_max');prior{ip}.cells_N_max=10;end
if ~isfield(prior{ip},'cells_N')
    if isfield(prior{ip},'cells_N_max')&&isfield(prior{ip},'cells_N_min')
        prior{ip}.cells_N=ceil(rand(1)*(prior{ip}.cells_N_max-prior{ip}.cells_N_min)+prior{ip}.cells_N_min);
    else
        prior{ip}.cells_N=10;
    end   
end
if prior{ip}.cells_N_min>prior{ip}.cells_N;  prior{ip}.cells_N_min=prior{ip}.cells_N; end
if prior{ip}.cells_N_max<prior{ip}.cells_N;  prior{ip}.cells_N_max=prior{ip}.cells_N; end


if ~isfield(prior{ip},'cells_use');
    prior{ip}.cells_use=randomsample(prior{ip}.cells_N_max,prior{ip}.cells_N);
end

%disp(sprintf('cells_N=%g',prior{ip}.cells_N));

if (nargin>1)&(prior{ip}.seq_gibbs.step>0)
    % perturb cells_N, only if STEP length > 0   
    step=prior{ip}.seq_gibbs.step;
    
    try
        step=max([1 ceil(prior{ip}.seq_gibbs.step)]);
    catch
        step=5;
    end
    direc=sign(rand(1)-.5);
    %disp(sprintf('%s: N=%d, moving=%d,',mfilename,prior{ip}.cells_N,step*direc));
    
    N_old=prior{ip}.cells_N;
    
    prior{ip}.cells_N=prior{ip}.cells_N+step*direc;
    if prior{ip}.cells_N>prior{ip}.cells_N_max; 
        prior{ip}.cells_N=N_old;
        direc=0;
    end
    if prior{ip}.cells_N<prior{ip}.cells_N_min; 
        prior{ip}.cells_N=N_old;
        direc=0;
    end
     
    
    % adjust cells_use
    % 
    if direc==-1
        % Removing cells
        % remove randomly selected cells
        ii=randomsample(prior{ip}.cells_use,step); 
        % update list of cells
        prior{ip}.cells_use=setxor(prior{ip}.cells_use,ii);
        
    elseif direc==1
        i_new=setxor(1:prior{1}.cells_N_max,prior{ip}.cells_use);
        if length(i_new)==1
            i_cell=i_new;
        else
            i_cell=randomsample(i_new,step);
        end
        prior{ip}.cells_use=[prior{ip}.cells_use,i_cell];
    elseif direc==0
        %% disp(sprintf('%s: staying, direc=%d',mfilename,direc))
    end
       
    
end


if ~(length(prior{ip}.cells_use)==prior{ip}.cells_N)
    % this means prior{ip}.cells_N must have been set by another 'prior'
    % variable, and then we mush adjust the number of used cells..
%    disp(sprintf('ncells_use=%d, N=%d',length(prior{ip}.cells_use),prior{ip}.cells_N));
    
    if length(prior{ip}.cells_use)>prior{ip}.cells_N
        if prior{ip}.cells_N==1,
            prior{ip}.cells_use=prior{ip}.cells_use.ceil(rand(1)*length(prior{ip}.cells_use));
        else
            prior{ip}.cells_use=randomsample(prior{ip}.cells_use,prior{ip}.cells_N);
        end
    else
        i_possible=setxor(1:prior{ip}.cells_N_max,prior{ip}.cells_use);
        N=prior{ip}.cells_N-length(prior{ip}.cells_use);
        cells_use_extra=randomsample(i_possible,N);
        prior{ip}.cells_use=[prior{ip}.cells_use,cells_use_extra];
    end
%   disp(sprintf('ncells_use=%d, N=%d',length(prior{ip}.cells_use),prior{ip}.cells_N));
end

% centers
if ~isfield(prior{ip},'cells_center');
    if prior{ip}.ndim>=1;
        xlim=[min(prior{ip}.x) max(prior{ip}.x)];
        prior{ip}.cells_center(:,1)=rand(prior{ip}.cells_N_max,1)*diff(xlim)+xlim(1);
    end
    if prior{ip}.ndim>=2;
        ylim=[min(prior{ip}.y) max(prior{ip}.y)];
        prior{ip}.cells_center(:,2)=rand(prior{ip}.cells_N_max,1)*diff(ylim)+ylim(1);
    end
    if prior{ip}.ndim>=3;
        zlim=[min(prior{ip}.z) max(prior{ip}.z)];
        prior{ip}.cells_center(:,3)=rand(prior{ip}.cells_N_max,1)*diff(zlim)+zlim(1);
    end
    
end
% cell value
if ~isfield(prior{ip},'cells_value');
    prior{ip}.cells_value=[1:prior{ip}.cells_N_max]';
end


%% compute nn map
i_use=prior{ip}.cells_use;
if prior{ip}.ndim==1;
    m_propose{ip} = interp1(prior{ip}.cells_center(i_use,1),prior{ip}.cells_value(i_use),prior{ip}.xx,'nearest','extrap');      
elseif prior{ip}.ndim==2;
    m_propose{ip} = griddata(prior{ip}.cells_center(i_use,1),prior{ip}.cells_center(i_use,2),prior{ip}.cells_value(i_use),prior{ip}.xx,prior{ip}.yy,'nearest');
elseif prior{ip}.ndim==3;
    m_propose{ip} = griddata(prior{ip}.cells_center(i_use,1),prior{ip}.cells_center(i_use,2),prior{ip}.cells_center(i_use,2),prior{ip}.cells_value(i_use),prior{ip}.xx,prior{ip}.yy,prior{ip}.yy,'nearest');
    % ??
end
    
if isempty(m_propose{ip});
    keyboard;
end

if nargin==1
    prior{ip}=rmfield(prior{ip},'cells_center');
    prior{ip}=rmfield(prior{ip},'cells_value');
    %prior{ip}=rmfield(prior{ip},'cells_use');
end



