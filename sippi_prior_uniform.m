% sippi_prior_uniform : Uniform prior for SIPPI
%
%
%% Example 1D uniform
%   ip=1;
%   prior{ip}.type='uniform';
%   prior{ip}.min=10;
%   prior{ip}.max=25;
%   [m,prior]=sippi_prior_uniform(prior);
%   sippi_plot_prior_sample(prior);
%
%% Example 10D uniform
%   ip=1;
%   prior{ip}.type='uniform';
%   prior{ip}.x=1:1:10; % As dimensions are uncorrelated, only the length
%                       % of prior{ip}.x matters, not its actual values.
%   prior{ip}.min=10;
%   prior{ip}.max=25;
%   [m,prior]=sippi_prior_uniform(prior);
%   sippi_plot_prior_sample(prior);
%
%% Sequential Gibbs sampling
%   prior{1}.seq_gibbs.step=.1;
%   for i=1:1000;
%       [m,prior]=sippi_prior(prior,m);
%       mm(i)=m{1};
%   end
%   subplot(1,2,1);plot(mm);
%   subplot(1,2,2);hist(mm);
%
% TMH/2014
%
% See also: sippi_prior_init, sippi_prior
%
function [m_propose,prior]=sippi_prior_uniform(prior,m_current,ip);

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

if nargin == 1
    if prior{1}.dim(3)>1
        prior{ip}.randn=randn([prior{ip}.dim(2) prior{ip}.dim(1) prior{ip}.dim(3)]);
    elseif prior{1}.dim(2)>1
        prior{ip}.randn=randn([prior{ip}.dim(2) prior{ip}.dim(1)]);
    else 
        prior{ip}.randn=randn(prior{ip}.dim);
    end
else
    % PERTURB
    
    if prior{ip}.seq_gibbs.type==1       
        % BOX TYPE RESIM
        new_randn=prior{ip}.randn;
        
        if (prior{1}.ndim)<=1
            step=prior{ip}.seq_gibbs.step;
            i_rand=randsample(size(new_randn,1),max([1 ceil(step)]));
            new_randn(i_rand)=randn(size(i_rand));            
        else
            disp(sprintf('%s: not yet implemeted seq gibbs type',mfilename))
        end
        
        %keyboard
        
        theta=.1*pi/2;
        prior{ip}.randn(i_rand) = prior{ip}.randn(i_rand)*cos(theta) + new_randn(i_rand)*sin(theta) ;
        %prior{ip}.randn=new_randn;
        
    elseif prior{ip}.seq_gibbs.type==2
        %
        
        % NEXT 7 LINES REPLACED BY 8TH LINE
        %if (prior{1}.ndim)<=1
        %    new_randn=randn(prior{ip}.dim);
        %elseif (prior{1}.ndim)==2
        %    new_randn=randn([prior{ip}.dim(2) prior{ip}.dim(1)]);
        %elseif (prior{1}.ndim)==3
        %    new_randn=randn([prior{ip}.dim(2) prior{ip}.dim(1) prior{ip}.dim(3)]);
        %end
        new_randn = randn(size(m_current{ip}));
        
        theta=90*(prior{ip}.seq_gibbs.step)*pi/180;
        prior{ip}.randn = prior{ip}.randn*cos(theta) + new_randn*sin(theta) ;
        
    end
end

m_propose{ip}=normcdf(prior{ip}.randn,0,1)*(prior{ip}.max-prior{ip}.min)+prior{ip}.min;
