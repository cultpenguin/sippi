% sippi_forward_linefit : line fit forward solition 
%
% [d,forward,prior,data]=sippi_forward_linefit(m,forward,prior,data);
%
function [d,forward,prior,data]=sippi_forward_linefit(m,forward,prior,data);

if length(m)==1;
    d{1}=forward.x*m{1};
else
    d{1}=forward.x*m{1}+m{2};
end
