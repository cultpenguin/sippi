% sippi_forward_fmm: return frequency distribution from a 1D/2D model
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
% Lange et al., 2012. A Frequency Matching Method: Solving Inverse Problems by Use of Geologically Realistic Prior Information. 
% Mathematical Geosciences, 44(7), 783-803, 2012. doi:10.1007/s11004-012-9417-2.
%
% See also: sippi_likelihood_fmm, frequency_matching
%

function [d,forward,prior,data]=sippi_forward_fmm(m,forward,prior,data,id,im)

if nargin<6, im=1;end
if nargin<5, id=1;end
if nargin<4, data{1}.null='';end

if nargin<2, forward.null='';end



if ~isfield(forward,'fmm');
    forward.fmm.null='';
end

if ~isfield(forward.fmm,'N');
    forward.fmm.N=2;
end

[d_hist,forward.fmm]=frequency_matching(round(m{1}),forward.fmm);
d{id}=d_hist(:);
%d{id}=d_hist(1:forward.fmm.N);
%d{id}=d_hist(:)./sum(d_hist(:));
