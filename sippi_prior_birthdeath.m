% sippi_prior_birthdeath
%
% Call: 
%     [m,prior]=sippi_prior_birthdeath(prior,m_current,im)
%
%
% prior{im}.type='birthdeath'
% prior{im}.N_layers_min % min number of layers
% prior{im}.N_layers_max % max number of layers
% prior{im}.v_min        % min value in layer
% prior{im}.v_max        % max value in layer
% prior{im}.z_interface_step % step in percentage of y-axis range, when
%                            % moving layer
%
% prior{im}.p_lev=[p_birth p_death p_move p_value]
%                 [1/6 1/6 1/6 1/2];
%                    p_birth, probability of birth of layer
%                    p_death, probability of death of layer
%                    p_move, probability of movement of layer boundary
%                    p_value, probability of perturbing value in each layer
%                    (step--> prior{im}.seq_gibbs.step)
%

function [m,prior]=sippi_prior_birthdeath(prior,m_current,im);

% if nargin == 01
%     for i=1:100
%         oimdisp(1)
%         [m,prior]=sippi_prior_birthdeath;
%         plot(m{1}+i*5,prior{1}.x);
%         hold on
%     end
%     hold off
%     return
% end

if nargin<3, im=1; end
if nargin<1, prior{im}.type='birthdeath';end

% force uncondtional
if nargin<2,
    if isfield(prior{im},'z_interface');
        prior{im}=rmfield(prior{im},'z_interface');
    end
    if isfield(prior{im},'v_interface');
        prior{im}=rmfield(prior{im},'v_interface');
    end
    if isfield(prior{im},'N_layers');
        prior{im}=rmfield(prior{im},'N_layers');
    end
end
    
%
if ~isfield(prior{im},'x');prior{im}.y=0:1:124;end

if ~isfield(prior{im},'N_layers_min');prior{im}.N_layers_min=1;end
if ~isfield(prior{im},'N_layers_max');prior{im}.N_layers_max=5;end
if ~isfield(prior{im},'N_layers');
    l_r = randi(prior{im}.N_layers_max-prior{im}.N_layers_min+1);
    prior{im}.N_layers=prior{im}.N_layers_min+l_r-1;
end

if ~isfield(prior{im},'v_min');prior{im}.v_min=-1;end
if ~isfield(prior{im},'v_max');prior{im}.v_max=3;end

% set z_interface, v_interface
if ~isfield(prior{im},'z_interface');
    y0=min(prior{im}.y);
    wy=max(prior{im}.y)-min(prior{im}.y);
    prior{im}.z_interface=sort(rand(1,prior{im}.N_layers-1)*wy+y0);
end
if ~isfield(prior{im},'z_interface_step');
    prior{im}.z_interface_step=0.02; % step length moving boundary
end



if ~isfield(prior{im},'v_interface');
    prior{im}.v_min=-1;
    prior{im}.v_max=3;
    wv=prior{im}.v_max-prior{im}.v_min;
    prior{im}.v_interface=rand(1,prior{im}.N_layers)*wv+prior{im}.v_min;
end

if ~isfield(prior,'init')
    prior=sippi_prior_init(prior);
end


%% birth deaths
if ~isfield(prior{im},'p_lev');
    prior{im}.p_lev=[1/6 1/6 1/6 1/2];
end
prior{im}.p_lev=prior{im}.p_lev./sum(prior{im}.p_lev);

pcum=cumsum(prior{im}.p_lev);
r=rand(1);
if r<pcum(1)&&(prior{im}.N_layers<prior{im}.N_layers_max);
    % birth
    sippi_verbose('birth',2)
    
    wv=prior{im}.v_max-prior{im}.v_min;
    v_interface_new = rand(1)*wv+prior{im}.v_min;
    
    y0=min(prior{im}.y);
    wy=max(prior{im}.y)-min(prior{im}.y);
    
    z_interface_new = rand(1)*wy+y0;

    sippi_verbose(sprintf('nz=%d',length(prior{im}.z_interface)),3)
    sippi_verbose(sprintf('nv=%d',length(prior{im}.v_interface)),3)

    interfaces_org = [[0,prior{im}.z_interface];[prior{im}.v_interface]];
    interfaces_new = [[0,prior{im}.z_interface, z_interface_new];[prior{im}.v_interface, v_interface_new]];
    interfaces_sort= sortrows(interfaces_new',1)';
    
    prior{im}.z_interface=interfaces_sort(1,2:end);
    prior{im}.v_interface=interfaces_sort(2,:);
    prior{im}.N_layers = prior{im}.N_layers + 1;
    prior{im}.N_layers = length(prior{im}.v_interface);
      
    sippi_verbose(sprintf('nz=%d',length(prior{im}.z_interface)),3)
    sippi_verbose(sprintf('nv=%d',length(prior{im}.v_interface)),3)
    
    
elseif r<pcum(2)&&(prior{im}.N_layers>prior{im}.N_layers_min);
    % death
    sippi_verbose('death',2)
    
    N_interfaces = prior{im}.N_layers-1;
    idel = randi(N_interfaces);
    sippi_verbose(sprintf('removing interface %d of %d',idel,N_interfaces),2)
    v_interface = prior{im}.v_interface(setxor(1:prior{im}.N_layers,idel+1));
    z_interface = prior{im}.z_interface(setxor(1:N_interfaces,idel));
    prior{im}.v_interface = v_interface;
    prior{im}.z_interface = z_interface;
    prior{im}.N_layers = length(prior{im}.v_interface);
    
    
elseif r<pcum(3);
    % move
    sippi_verbose('move',2)
    
    N_interfaces = prior{im}.N_layers-1;
    imove = randi(N_interfaces);
    
    move_step = prior{im}.z_interface_step;
    dz = max(prior{1}.y)-min(prior{1}.y);
    
    z_interface = prior{im}.z_interface;   
    z_interface(imove) = z_interface(imove) + randn(1)*move_step*dz;
    
    if (z_interface(imove)>min(prior{1}.y))&&(z_interface(imove)<max(prior{1}.y))
        prior{im}.z_interface = z_interface;
    end
    
    
else
    % resistvity
    sippi_verbose('resistivity',2)
    prior{im}.v_interface = prior{im}.v_interface + randn(size(prior{im}.v_interface)).*prior{im}.seq_gibbs.step;
    imax=find(prior{im}.v_interface>prior{im}.v_max);
    imin=find(prior{im}.v_interface<prior{im}.v_min);
    prior{im}.v_interface(imin)=prior{im}.v_min;
    prior{im}.v_interface(imax)=prior{im}.v_max;
end
    

sippi_verbose(sprintf('Nl=%d',prior{im}.N_layers),2)
    

%% remove thin layers


%% build model
m{im}=ones(length(prior{im}.y),1).*prior{im}.v_interface(1);
for i=1:(prior{im}.N_layers-1);
    ii=find(prior{im}.y>prior{im}.z_interface(i));
    m{im}(ii)=prior{im}.v_interface(i+1);
end



%% return proposed model



    
    
    