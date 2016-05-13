% sippi_forward_fmm: reutrn frequency distribution
%
% Call :
%  [d,forward,prior,data]=sippi_forward_fmm(m,forward,prior,data,id,im)
%
% See also frequency_matching
%
%  ip=1;
%     prior{ip}.type='mps';
%     prior{ip}.method='mps_snesim';
%     prior{ip}.x=1:1:80;
%     prior{ip}.y=1:1:80;
%     prior{ip}.ti=channels;
%     m=sippi_prior(prior);
%
%     [d,forward]=sippi_forward_fmm(m)
%
%     forward.forward_function='sippi_forward_fmm';
%     [d,forward]=sippi_forward(m,forward)
% 
%
function [d,forward,prior,data]=sippi_forward_fmm(m,forward,prior,data,id,im)

if nargin<6, im=1;end
if nargin<5, id=1;end
if nargin<4, data{1}.null='';end

if nargin<2, forward.null='';end

if ~isfield(forward,'fmm');
    forward.fmm.null='';
end

[d_hist,forward.fmm]=frequency_matching(m{1},forward.fmm);
d{id}=d_hist(:);
