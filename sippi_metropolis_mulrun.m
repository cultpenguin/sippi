% sippi_metropolis_mulrun: multiple (independent) Metropolis chains in parallel
%
% Runs multiple independent Metropolis chains.
% If the Matlab parallel toolbox is available, each chain will be run 
% on a different thread. 
% This should provide close to linear speedup with the number of avilable
% threads.
% 
% To start 8 independent Metropolis samplers, use:
%   options.nruns=8;
%   [options_mul]=sippi_metropolis_mulrun(data,prior,forward,options);
%
%
% See also: sippi_metropolis
%
function [options_mul]=sippi_metropolis_mulrun(data,prior,forward,options)
if nargin<3, 
    help(mfilename);
    options_mul='';
    return
end

if nargin<4, options.nruns=2;end

if ~isfield(options,'nruns')
    options.nruns=2;
end

parfor i=1:options.nruns
    o=options;
    
    if isfield(options,'txt');
        o.txt=sprintf('%s_r%02d',options.txt,i);
    else
        o.txt=sprintf('r%02d',i);
    end
    
    % sippi_metropolis Extended Metropolis sampling in SIPPI
    [options_mul{i}]=sippi_metropolis(data,prior,forward,o);
end
        