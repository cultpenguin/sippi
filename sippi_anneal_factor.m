% sippi_anneal_factor : compute simple noise multiplication factor for
% annealing type sampling
%
% See also sippi_metropolis, sippi_anneal_adjust_noise
%
function [fac,mcmc]=sippi_anneal_factor(mcmc,i,prior);

if nargin==0;
    mcmc.anneal.i_end=20000;
    i=1:1:2*mcmc.anneal.i_end;
    [fac,mcmc]=sippi_anneal_factor(mcmc,i);
    plot(i,fac,'k*')
    xlabel('Iteration number')
    ylabel('Noise amplification')
    return
end

if ~isfield(mcmc,'anneal');mcmc.anneal.null='';end
if ~isfield(mcmc.anneal,'type');mcmc.anneal.type='exp';end
%if ~isfield(mcmc.anneal,'type');mcmc.anneal.type='linear';end
if ~isfield(mcmc.anneal,'fac_begin');mcmc.anneal.fac_begin=20;end
if ~isfield(mcmc.anneal,'fac_end');mcmc.anneal.fac_end=1;end
if ~isfield(mcmc.anneal,'i_begin');mcmc.anneal.i_begin=1;end
if ~isfield(mcmc.anneal,'i_end');
    try
        mcmc.anneal.i_end=1e+9;
        for ip=1:length(prior);
            mcmc.anneal.i_end=min([prior{ip}.seq_gibbs.i_update_step_max.*.5 mcmc.anneal.i_end]);
        end
    catch
        mcmc.anneal.i_end=1;
    end
    
end

if strcmp(mcmc.anneal.type,'linear');
    % linear
    if i>mcmc.anneal.i_end;
        fac=mcmc.anneal.fac_end;
    elseif i<mcmc.anneal.i_begin;
        fac=mcmc.anneal.fac_begin;
    else
        delta_i=(i-mcmc.anneal.i_begin)./(mcmc.anneal.i_end-mcmc.anneal.i_begin);
        delta_fac=(mcmc.anneal.fac_end-mcmc.anneal.fac_begin);
        fac_change=delta_i*delta_fac;
        fac = mcmc.anneal.fac_begin+fac_change;
    end
    
    if length(i)>1;
        i1=find(i>mcmc.anneal.i_end);
        fac(i1)=mcmc.anneal.fac_end;
        
        i2=find(i<mcmc.anneal.i_begin);
        fac(i2)=mcmc.anneal.fac_begin;
    end
    
elseif strcmp(mcmc.anneal.type,'exp');
    %%
    x1=mcmc.anneal.i_begin;
    x2=mcmc.anneal.i_end;
    y1=mcmc.anneal.fac_begin;
    y2=mcmc.anneal.fac_end;
    
    b=(y1/y2)^(1/(x1-x2));
    a=y1/b.^(x1);
    
    mcmc.anneal.a=a;
    mcmc.anneal.b=b;
    
    fac = mcmc.anneal.a.*mcmc.anneal.b.^(i);
    
    i1=find(i>mcmc.anneal.i_end);
    fac(i1)=mcmc.anneal.fac_end;
    
    i2=find(i<mcmc.anneal.i_begin);
    fac(i2)=mcmc.anneal.fac_begin;
    
    
end

