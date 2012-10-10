% sippi_adjust_step_size Adjust step length length for Metropolis sampler in SIPPI
%  
% Call : 
%   step=sippi_adjust_step_size(step,P_average,P_target);
%
% step : current step 
% P_current : Current acceptance ratio
% P_target  : preferred acceptance ratio (def=0.3);
%
% See also sippi_compute_acceptance_rate, sippi_prior_set_steplength
%
function step=sippi_adjust_step_size(step,P_current,P_target);

if nargin<3
    P_target=0.3;
end

step = step.*(P_current/P_target).^(.25);
