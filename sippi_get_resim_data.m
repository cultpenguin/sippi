% sippi_get_resim_data: Get conditional data for resimulation
%
% d_cond=sippi_get_resim_data(m_current,prior,ip);
% 
% c_cond [n_cond,4]: col1: x, col2: y, col4: z, col4: d
%
% See also sippi_prior
%
function d_cond=sippi_get_resim_data(m_current,prior,ip);
if nargin<3
    ip=1;
end

step=prior{ip}.seq_gibbs.step;
if prior{ip}.seq_gibbs.type==1;
    if ~isfield(prior{ip}.seq_gibbs,'n_ite');
        prior{ip}.seq_gibbs.n_ite=1;
    end
    n_ite=prior{ip}.seq_gibbs.n_ite;
        
    if length(step)<2, step(2:3)=step(1);end
    if length(step)<3, step(3)=step(2);end
  
    ih=[];
    used=prior{ip}.xx.*0+1;
    for j=1:n_ite
        % 
        pos(1)=min(prior{ip}.x)+rand(1)*(max(prior{ip}.x)-min(prior{ip}.x));
        pos(2)=min(prior{ip}.y)+rand(1)*(max(prior{ip}.y)-min(prior{ip}.y));
        pos(3)=min(prior{ip}.z)+rand(1)*(max(prior{ip}.z)-min(prior{ip}.z));

        used(find( (abs(prior{ip}.xx-pos(1))<step(1)) & (abs(prior{ip}.yy-pos(2))<step(2)) ))=0;
    end
    ih=find(used);
    %if j>1; ih=unique(ih);end
    d_cond=[prior{ip}.xx(ih(:)) prior{ip}.yy(ih(:)) prior{ip}.zz(ih(:)) m_current{ip}(ih(:))];

    
elseif prior{ip}.seq_gibbs.type==2
    % RANDOM SELECTION OF MODEL PARAMETERS FOR RESIM
    N=prod(size(m_current{ip}));
    n_resim=prior{1}.seq_gibbs.step(1);
    if n_resim<=1;
        % if n_resim is less than one
        % n_resim defines the fraction of N to use
        n_resim=ceil(n_resim.*N);
    end
    n_cond=N-n_resim;
    ih=randomsample(N,n_cond);
    d_cond=[prior{ip}.xx(ih(:)) prior{ip}.yy(ih(:)) prior{ip}.zz(ih(:)) m_current{ip}(ih(:))];
end
