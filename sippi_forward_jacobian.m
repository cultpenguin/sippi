% sippi_forward_jacobian: Compute jacobian / partial derivative
%
% Call: 
%   J=sippi_forward_jacobian(m,forward,prior);
%   J=sippi_forward_jacobian(m,forward,prior,dm,);
%   J=sippi_forward_jacobian(m,forward,prior,dm,used);
%
% In:
%   m: SIPPI model, as in m=sippi_prior(prior(;
%   forward: SIPPI forward structure
%   prior: SIPPI prior structure
%   dm: pertubation to model parameter (def: dm=0.001*mean(m{1}(:)))
%   used: Compute J, for data number 'used' (def: used=1)
% Out:
%   J [ny,nx,nz]: Jacobian matrix of partial derivatives
%                (d(m_i)-d(m_i+dm))/dm;
%
%
function J=sippi_forward_jacobian(m,forward,prior,dm,used);

id=1;
ip=1;

if ~isfield(prior,'init')
    prior=sippi_prior_init(prior);
end

J=zeros([ prior{ip}.dim(2) prior{ip}.dim(1) prior{ip}.dim(3) ]);

if isempty(m);
    m=sippi_prior(prior);
end

d=sippi_forward(m,forward,prior);

if nargin<4
    dm=0.001*mean(m{1}(:));
end
    

if nargin<5
    used=1;
end


for ix=1:prior{ip}.dim(1);
    progress_txt([ix],[prior{ip}.dim(1)])
    for iy=1:prior{ip}.dim(2);
        %progress_txt([ix iy],[prior{ip}.dim(1),prior{ip}.dim(2)])
        for iz=1:prior{1}.dim(3);
            
            m_p=m;
            m_p{1}(iy,ix,iz)=m_p{ip}(iy,ix,iz)+dm;
            d_p=sippi_forward(m_p,forward,prior);
            
            J(iy,ix,iz)=(d{1}(used)-d_p{id}(used))./dm;
            
        end
        
    end
end


