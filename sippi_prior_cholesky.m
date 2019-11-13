% sippi_prior_cholesky : Cholesky type Gaussian prior for SIPPI
%
%% Example:
%   ip=1;
%   prior{ip}.type='cholesky';
%   prior{ip}.m0=10;
%   prior{ip}.Cm='.001 Nug(0) + 1 Gau(10)';
%   prior{ip}.x=0:1:100;linspace(0,100,20);
%   prior{ip}.y=0:1:50;linspace(0,33,30);
%   [m,prior]=sippi_prior_cholesky(prior);
%   sippi_plot_prior(prior,m);
%
%% Sequential Gibbs sampling
%   prior{1}.seq_gibbs.step=.1;
%   for i=1:100;
%       [m,prior]=sippi_prior_cholesky(prior,m);
%       sippi_plot_prior(prior,m);
%       caxis([8 12]);drawnow;
%   end
%
%% Prior covariance model
% The prior covarince model can be setup using
%   prior{ip}.m0=10;
%   prior{ip}.Cm='.001 Nug(0) + 1 Gau(10)';
% or
%   prior{ip}.m0=10;
%  and the 'Cmat' variable 'prior{ip}.Cmat' which much the contain a full
%  nd X nd size covariance matrix.
%  (it is computed the first the sippi_prior_cholesky is called)
%
% See also: gaussian_simulation_cholesky, sippi_sequential_gibbs_resim
%
function [m_propose,prior]=sippi_prior_cholesky(prior,m_current,ip);

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

if (isfield(prior{ip},'d_target'))&(~isfield(prior{ip},'o_nscore'))
    % UPDATE PRIOR STRUCTURE TO USE TARGET DISTRIBUTION
    d_min=min(prior{ip}.d_target);
    d_max=max(prior{ip}.d_target);
    [d_nscore,o_nscore]=nscore(prior{ip}.d_target,1,1,d_min,d_max,0);
    prior{ip}.o_nscore=o_nscore;
end

%% Remove Rand Number if they are not needed
if (nargin==1)&isfield(prior{ip},'z_rand');
    prior{ip}=rmfield(prior{ip},'z_rand');
end

%% Sequential Gibbs?
if nargin>1
    
    z_cur=prior{ip}.z_rand;
    z_init=z_cur;
    z_new=randn(size(z_cur));
    
    i_resim = sippi_sequential_gibbs_resim(prior,ip);
    
    z_cur(i_resim)=randn(size(z_cur(i_resim)));    
    %% linear combinartion of the perturbed paramaters
    if (isfield(prior{ip}.seq_gibbs,'gradual'))
        if prior{ip}.seq_gibbs.gradual<1
            
            if exist('gaussian_linear_combine','file')
                i_perturbed=i_resim;%find((z_init-z_cur)~=0);
                z_cur(i_perturbed) = gaussian_linear_combine(z_init(i_perturbed),z_cur(i_perturbed),prior{ip}.seq_gibbs.gradual,0);
            end
        end
    end

    prior{ip}.z_rand=z_cur;

end



%%
nsim=1;
if isfield(prior{ip},'L')
    is_chol=1;
    if isfield(prior{ip},'z_rand')
        [z,D,prior{ip}.L,prior{ip}.z_rand]=gaussian_simulation_cholesky(0,prior{ip}.L,nsim,is_chol,prior{ip}.z_rand);
    else
        [z,D,prior{ip}.L,prior{ip}.z_rand]=gaussian_simulation_cholesky(0,prior{ip}.L,nsim,is_chol);
    end
else
    % Compute Cov matrix
    if ~isfield(prior{ip},'Cmat')
        if prior{ip}.ndim==1;
            prior{ip}.Cmat=precal_cov([prior{ip}.xx(:)],[prior{ip}.xx(:)],prior{ip}.Cm);
        elseif prior{ip}.ndim==2;
            prior{ip}.Cmat=precal_cov([prior{ip}.xx(:) prior{ip}.yy(:)],[prior{ip}.xx(:) prior{ip}.yy(:)],prior{ip}.Cm);
        elseif prior{ip}.ndim==3;
            prior{ip}.Cmat=precal_cov([prior{ip}.xx(:) prior{ip}.yy(:) prior{ip}.zz(:)],[prior{ip}.xx(:) prior{ip}.yy(:) prior{ip}.zz(:)],prior{ip}.Cm);
        end        
    end

    is_chol=0;
    if isfield(prior{ip},'z_rand')
        [z,D,prior{ip}.L,prior{ip}.z_rand]=gaussian_simulation_cholesky(0,prior{ip}.Cmat,nsim,is_chol,prior{ip}.z_rand);
    else
        [z,D,prior{ip}.L,prior{ip}.z_rand]=gaussian_simulation_cholesky(0,prior{ip}.Cmat,nsim,is_chol);
    end
end

% return data to m_propose
if prior{ip}.ndim==1;
    D=D';
elseif prior{ip}.ndim==2;
    D=reshape(z,prior{ip}.dim(2),prior{ip}.dim(1));
elseif prior{ip}.ndim==3;
    D=reshape(z,prior{ip}.dim(2),prior{ip}.dim(1),prior{ip}.dim(3));
end
m_propose{ip}=D;
if isfield(prior{ip},'o_nscore')

    if ~isstruct(prior{ip}.Va);
        prior{ip}.Va=deformat_variogram(prior{ip}.Va);
    end
    Va_par=prior{ip}.Va;
    gvar=sum([Va_par.par1]);
    D=D./sqrt(gvar);

    m_propose{ip}=inscore(D,prior{ip}.o_nscore)+prior{ip}.m0;
else    
    m_propose{ip}=D+prior{ip}.m0;
end
