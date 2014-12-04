% sippi_prior_voronoi: 
%
% TMH/2014
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
if ~isfield(prior{ip},'cells_N');
    prior{ip}.cells_N=10;
end

%if nargin<2
%    if isfield(prior{ip},'cells_center');
%        prior{ip}=rmfield(prior{ip},'cells_center');
%    end
%    if isfield(prior{ip},'cells_value');
%        prior{ip}=rmfield(prior{ip},'cells_value');
%    end
%end

% centers
if ~isfield(prior{ip},'cells_center');
    if prior{ip}.ndim>=1;
        xlim=[1 100];
        prior{ip}.cells_center(:,1)=rand(prior{ip}.cells_N,1)*diff(xlim)+xlim(1);
    end
    if prior{ip}.ndim>=2;
        ylim=[1 100];
        prior{ip}.cells_center(:,2)=rand(prior{ip}.cells_N,1)*diff(ylim)+ylim(1);
    end
    if prior{ip}.ndim>=2;
        zlim=[1 100];
        prior{ip}.cells_center(:,2)=rand(prior{ip}.cells_N,1)*diff(zlim)+zlim(1);
    end
    
end

% cell value
if ~isfield(prior{ip},'cells_value');
    prior{ip}.cells_value=[1:prior{ip}.cells_N]';
end


%% compute nn map
if prior{ip}.ndim==1;
    
elseif prior{ip}.ndim==2;
    m_propose{ip} = griddata(prior{ip}.cells_center(:,1),prior{ip}.cells_center(:,2),prior{ip}.cells_value,prior{ip}.xx,prior{ip}.yy,'nearest');
elseif prior{ip}.ndim==3;
    
end
    



