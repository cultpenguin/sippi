% sippi_forward : forward wrapper
%
% Call : 
%  [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id,im)
%
%
function [d,forward,prior,data]=sippi_forward(m,forward,prior,data,id,im)

if nargin<5, id=1;end
if nargin<6, im=1;end

% TRAVEL TIME FORWARD
[d,forward,prior,data]=sippi_forward_traveltime(m,forward,prior,data,id,im);


