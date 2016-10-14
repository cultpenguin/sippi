% sippi_forward_linefit Line fit forward solver for SIPPI 
%
% [d,forward,prior,data]=sippi_forward_linefit(m,forward,prior,data);
% [d,forward,prior,data]=sippi_forward_linefit(m,forward);
%
function [d,forward,prior,data]=sippi_forward_linefit(m,forward,prior,data);

if length(m)==1;
    d{1}=forward.x.*0 + m{1};
elseif length(m)==2;
    d{1}=forward.x*m{2}+m{1};
else
    d{1}=forward.x.^2*m{3}+forward.x*m{2}+m{1};
end
