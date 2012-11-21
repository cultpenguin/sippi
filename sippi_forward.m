% sippi_forward Simple forward wrapper for SIPPI
%
% Assumes that the actual forward solver has been defined by
% forward.forward_function
%
% Call:
%   [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id,im)
%
function [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id,im)

%if nargin<4;    forward.null='';end
if nargin<4;    data{1}.null='';end
%if nargin<5;    id=1;end
%if nargin<6;    im=1;end


if isfield(forward,'forward_function');
    
    if nargin==2;
        [d,forward]=feval(forward.forward_function,m,forward);
    elseif nargin==3
        [d,forward,prior]=feval(forward.forward_function,m,forward,prior);
    elseif nargin==4
        [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior,data);
    elseif nargin==5
        [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior,data,id);
    else
        [d,forward,prior,data]=feval(forward.forward_function,m,forward,prior,data,id,im);
    end
else
    disp(sprintf('%s : No forward_function specified in ''forward'' structure',mfilename))
    d=[];
    
end

