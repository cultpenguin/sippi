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
% To manually set the 'local' cluster profile to allow using 4 threads
% use e.g.:
%    myCluster = parcluster('local');
%    myCluster.NumWorkers = 4;  % 'Modified' property now TRUE
%    saveProfile(myCluster);    % 'local' profile now updated,
%                               % 'Modified' property now FALSE   
%
% See also: sippi_metropolis
%
function [options_mul,data,prior,forward]=sippi_metropolis_mulrun(data,prior,forward,options)
if nargin<3, 
    help(mfilename);
    options_mul='';
    return
end

if nargin<4, options.nruns=2;end

if ~isfield(options,'nruns')
    options.nruns=2;
    sippi_verbose(sprintf('%s: setting options.nruns=%d', mfilename,options.nruns),-1);
end

sippi_verbose(sprintf('%s: Runnning sippi_metropolis %d times (as independent chains)', mfilename,options.nruns),-1);
parfor i=1:options.nruns
    o=options;
    
    if isfield(options,'txt');
        o.txt=sprintf('%s_r%02d',options.txt,i);
    else
        o.txt=sprintf('r%02d',i);
    end
    o.nruns=1; % Needed to avoid recursive call to sippi_metropolis_mulrun
    % sippi_metropolis Extended Metropolis sampling in SIPPI
    
    [options_mul{i}]=sippi_metropolis(data,prior,forward,o);
end
