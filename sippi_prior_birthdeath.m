function [m,prior]=sippi_prior_birthdeath(prior,m_current,im);

% if nargin == 01
%     for i=1:100
%         oimdisp(1)
%         [m,prior]=sippi_prior_voronoi;
%         plot(m{1}+i*5,prior{1}.x);
%         hold on
%     end
%     hold off
%     return
% end

if nargin<3, im=1; end
if nargin<1, prior{im}.type='birthdeath';end

% force uncondtional
if isfield(prior{im},'z_interface');
    prior{im}=rmfield(prior{im},'z_interface');
end
if isfield(prior{im},'v_interface');
    prior{im}=rmfield(prior{im},'v_interface');
end
if isfield(prior{im},'N_layers');
    prior{im}=rmfield(prior{im},'N_layers');
end

%
if ~isfield(prior{im},'x');prior{im}.y=0:1:124;end

if ~isfield(prior{im},'N_layers_min');prior{im}.N_layers_min=1;end
if ~isfield(prior{im},'N_layers_max');prior{im}.N_layers_max=5;end
if ~isfield(prior{im},'N_layers');
    l_r = randi(prior{im}.N_layers_max-prior{im}.N_layers_min+1);
    prior{im}.N_layers=prior{im}.N_layers_min+l_r-1;
end

% force uncondtional
if isfield(prior{im},'z_interface');
    prior{im}=rmfield(prior{im},'z_interface');
end
if isfield(prior{im},'v_interface');
    prior{im}=rmfield(prior{im},'v_interface');
end

% set z_interface, v_interface
if ~isfield(prior{im},'z_interface');
    x0=min(prior{im}.x);
    wx=max(prior{im}.x)-min(prior{im}.x);
    prior{im}.z_interface=sort(rand(1,prior{im}.N_layers-1)*wx+x0);
end

if ~isfield(prior{im},'val_interface');
    v_min=-1;
    v_max=3;
    wv=v_max-v_min;
    prior{im}.v_interface=rand(1,prior{im}.N_layers)*wv+x0;
end

%% birth deaths
if ~isfield(prior{im},'p_lev');
    prior{im}.p_lev=[1/6 1/6 1/6 1/2];
end
prior{im}.p_lev=prior{im}.p_lev./sum(prior{im}.p_lev);

pcum=cumsum(prior{im}.p_lev);
r=rand(1);
if r<pcum(1);
    % birth
elseif r<pcum(2);
    % death
elseif r<pcum(3);
    % move
else
    % nothing
end
    
% resistivity



%% remove thin layers


%% build model
m{im}=ones(length(prior{im}.y),1).*prior{im}.v_interface(1);
for i=1:(prior{im}.N_layers-1);
    ii=find(prior{im}.x>prior{im}.z_interface(i));
    m{im}(ii)=prior{im}.v_interface(i);
end



%% return proposed model



    
    
    