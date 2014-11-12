% sippi_forward Simple forward wrapper for SIPPI
%
% Assumes that the actual forward solver has been defined by
% forward.forward_function
%
% Call:
%   [d,forward,prior,data]=sippi_forward(m,forward)
%
% Optional: 
%   [d,forward,prior,data]=sippi_forward(m,forward,prior)
%   [d,forward,prior,data]=sippi_forward(m,forward,prior,data)
%   [d,forward,prior,data]=sippi_forward(m,forward,prior,data,options)
%
function [d,forward,prior,data,options]=sippi_forward(m,forward,prior,data,options)

if nargin<4;    data{1}.null='';end

% make sure to initilize the prior if it has not allready been done
% TMH: can this be ignored?
for ip=1:length(prior)
    if ~isfield(prior{ip},'init');
        prior=sippi_prior_init(prior);
    end
end

if isfield(forward,'forward_function');
    
    if nargin==2;
        [d,forward]=feval(forward.forward_function,m,forward);
    elseif nargin==3
        [d,forward,prior]=feval(forward.forward_function,m,forward,prior);
    elseif nargin==4
        [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior,data);
    elseif nargin==5
        [d,forward,prior,data,options]=feval(forward.forward_function,m,forward,prior,data,options);
    end
else
    disp(sprintf('%s : No forward_function specified in ''forward'' structure',mfilename))
    d=[];
    
end

