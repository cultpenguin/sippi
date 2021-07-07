% sippi_prior_layered
%
% Call:
%     [m,prior]=sippi_prior_layered(prior,m_current,im)
%
%
% Example 
%         im=1;
%         prior{im}.type = 'layered';
%         prior{im}.name = 'GrevieLayer';
%         prior{im}.x=z;
% 
%         j=0;
%         j=j+1;
%         prior{im}.p_rho{j}.name = 'PreCat';
%         prior{im}.p_rho{j}.type = 'uniform';
%         prior{im}.p_rho{j}.min = 20;prior{im}.p_rho{j}.max = 300;
%         prior{im}.p_thick{j}.name = sprintf('THICK_%02d',j);
%         prior{im}.p_thick{j}.type = 'uniform';
%         prior{im}.p_thick{j}.min = 20;
%         prior{im}.p_thick{j}.max = 40;
% 
%         j=j+1;
%         prior{im}.p_rho{j}.name = 'CatCat';
%         prior{im}.p_rho{j}.type = 'uniform';
%         prior{im}.p_rho{j}.min = 80;prior{im}.p_rho{j}.max = 120;
%         prior{im}.p_thick{j}.name = sprintf('THICK_%02d',j);
%         prior{im}.p_thick{j}.type = 'uniform';
%         prior{im}.p_thick{j}.min = 10;
%         prior{im}.p_thick{j}.max = 12;
% 
%         j=j+1;
%         prior{im}.p_rho{j}.name = 'Chalk';
%         prior{im}.p_rho{j}.type = 'uniform';
%         prior{im}.p_rho{j}.min = 1;prior{im}.p_rho{j}.max = 25;
%         prior{im}.p_thick{j}.name = sprintf('THICK_%02d',j);
%         prior{im}.p_thick{j}.type = 'uniform';
%         prior{im}.p_thick{j}.min = 0;
%         prior{im}.p_thick{j}.max = 300;
% 
%         j=j+1;
%         prior{im}.p_rho{j}.name = 'BASE';
%         prior{im}.p_rho{j}.type = 'uniform';
%         prior{im}.p_rho{j}.min = 40;
%         prior{im}.p_rho{j}.max = 40;
%
%
function [m,prior]=sippi_prior_layered(prior,m_current,im);

if nargin<3, im=1; end
if ~isfield(prior{1},'init')
    prior=sippi_prior_init(prior);
end
if ~isfield(prior{1},'name');prior{im}.name='rho';end
if ~isfield(prior{1},'is_cond');prior{im}.is_cond=0;end
if ~isfield(prior{1},'x')
    prior{1}.x = [0:1:300];
end

if ~isfield(prior{im},'cmap')
    % Colormap from Anne-sophie!
    prior{im}.cmap = [0,0,0;
        255,179,179;
        255,196,196;
        255,128,128;
        215,149,83;
        255,70,70;
        179,89,0;
        249,0,0;
        113,56,0;
        193,0,0;
        0,128,255;
        215,255,255;
        0,128,192;
        0,0,255,;
        0,128,0]/255;
end

if ~isfield(prior{1},'n_layers');prior{1}.n_layers=length(prior{1}.p_rho);end



if nargin>1
    try
    % update step lengths
    for i=1:length(prior{1}.p_rho)
        prior{1}.p_rho{i}.seq_gibbs.step = prior{1}.seq_gibbs.step;
    end
    for i=1:length(prior{1}.p_thick)
        if length(prior)>1
            prior{1}.p_thick{i}.seq_gibbs.step = prior{2}.seq_gibbs.step;
        else
            prior{1}.p_thick{i}.seq_gibbs.step = prior{1}.seq_gibbs.step;
        end
    end
    catch
        keyboard
    end
    
    [prior{1}.m_rho,prior{1}.p_rho] = sippi_prior(prior{1}.p_rho,prior{1}.m_rho);
    [prior{1}.m_thick,prior{1}.p_thick] = sippi_prior(prior{1}.p_thick,prior{1}.m_thick);
else
    [prior{1}.m_rho,prior{1}.p_rho] = sippi_prior(prior{1}.p_rho);
    [prior{1}.m_thick,prior{1}.p_thick] = sippi_prior(prior{1}.p_thick);
end
for i = 1:length(prior{1}.m_thick)
    m_t(i)=prior{1}.m_thick{i};
end
m_z = cumsum(m_t);

rho = ones(prior{1}.dim(1),1).*prior{1}.m_rho{end};
lith = ones(prior{1}.dim(1),1)*length(prior{1}.m_rho);


for i=(prior{1}.n_layers-1):-1:1
    ix=find(prior{im}.x<m_z(i));
    lith(ix)=i;
    rho(ix)=prior{1}.m_rho{i};
end

prior{im}.lith=lith;
prior{im}.thick=m_z;
if prior{im}.is_cond==1
    m{im}=log10(rho);
else
    m{im}=log10(1./rho);
end


