% sippi_anneal_temperature : compute annealing temperature for
% annealing type sampling
%
%   %% ANNEALING (TEMPERATURE AS A FUNTION OF ITERAITON NUMBER)
%    i % iteration number
%
%    mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
%    mcmc.anneal.i_end=100000; %  iteration number when annealing stops
%    mcmc.anneal.T_begin=5; % Start temperature for annealing 
%    mcmc.anneal.T_end=1; % End temperature for anneaing 
%
%    mcmc.anneal.type='exp';     % Exponential temperature change
%    mcmc.anneal.type='linear';  % Linear temperature change
% 
% Call
%   [T,mcmc]=sippi_anneal_temperature(i,mcmc);
%
% See also sippi_metropolis
%
function [T,mcmc]=sippi_anneal_temperature(i,mcmc,prior);

if nargin<2
    mcmc.anneal.null='';
end
if nargin==0;
    mcmc.anneal.i_end=20000;
    i=1:1:2*mcmc.anneal.i_end;
    [fac,mcmc]=sippi_anneal_temperature(mcmc,i);
    plot(i,fac,'k*')
    xlabel('Iteration number')
    ylabel('Noise amplification')
    return
end

if ~isfield(mcmc,'anneal');mcmc.anneal.null='';end

if isfield(mcmc.anneal,'fac_begin');mcmc.anneal.T_begin=mcmc.anneal.fac_begin;end
if isfield(mcmc.anneal,'fac_end');mcmc.anneal.T_begin=mcmc.anneal.fac_end;end

if ~isfield(mcmc.anneal,'type');mcmc.anneal.type='exp';end
%if ~isfield(mcmc.anneal,'type');mcmc.anneal.type='linear';end

if ~isfield(mcmc.anneal,'T_begin');mcmc.anneal.T_begin=5;end
if ~isfield(mcmc.anneal,'T_end');mcmc.anneal.T_end=1;end
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
    keyboard
    % linear
    if i>mcmc.anneal.i_end;
        T=mcmc.anneal.T_end;
    elseif i<mcmc.anneal.i_begin;
        T=mcmc.anneal.T_begin;
    else
        delta_i=(i-mcmc.anneal.i_begin)./(mcmc.anneal.i_end-mcmc.anneal.i_begin);
        delta_fac=(mcmc.anneal.T_end-mcmc.anneal.T_begin);
        T_change=delta_i*delta_fac;
        T = mcmc.anneal.T_begin+T_change;
    end
    
    if length(i)>1;
    
        i0=find(x1==x2);
        T(i0)=mcmc.anneal.T_end;
    
        i1=find(i>mcmc.anneal.i_end);
        T(i1)=mcmc.anneal.T_end;
        
        i2=find(i<mcmc.anneal.i_begin);
        T(i2)=mcmc.anneal.T_begin;
    end
    
elseif strcmp(mcmc.anneal.type,'cos');
    (cos(0:.001:pi)+1)*(mcmc.anneal.T_end-mcmc.anneal.T_begin)+mcmc.anneal.T_begin;
elseif strcmp(mcmc.anneal.type,'exp');
    %%
    x1=mcmc.anneal.i_begin;
    x2=mcmc.anneal.i_end;
    y1=mcmc.anneal.T_begin;
    y2=mcmc.anneal.T_end;
    
    b=(y1/y2)^(1/(x1-x2));
    a=y1/b.^(x1);
    
    mcmc.anneal.a=a;
    mcmc.anneal.b=b;
    
    T = mcmc.anneal.a.*mcmc.anneal.b.^(i);
    
    i0=find(x1==x2);
    T(i0)=mcmc.anneal.T_end;
    
    i1=find(i>mcmc.anneal.i_end);
    T(i1)=mcmc.anneal.T_end;
    
    i2=find(i<mcmc.anneal.i_begin);
    T(i2)=mcmc.anneal.T_begin;
    
end

