% sippi_forward_linear: 
%
% % options:
% forward.G : Linear forward operator. such that d[
%             d{id}=forward.G*m{im}(:)
%             if not set, forward.G=eye(prod(size(m{im})));
%
% forward.force_sparse [0]: Use forward.G as is (default)
%                     [1]: force forward.G to be treated as a sparse matrix
%
% %% Examples
% % Define an example of a prior model
%    forward.forward_function='sippi_forward_linear';
%    im=1;
%    prior{im}.type='FFTMA';
%    prior{im}.x=[0:1:100];
%    prior{im}.m0=0;
%    prior{im}.Va='1 Sph(10)';
% % Define the use of the linear forward solver
%    forward.forward_function='sippi_forward_linear';
%
% % Example 1: Identity, 
% %             by default an identity operator is used if the linear
% %             is not set
%   m=sippi_prior(prior);
%   tic;
%   d=sippi_forward(m,forward);
%   t1=toc;
%   figure(1);
%   subplot(1,2,1);plot(m{1});title('m')
%   subplot(1,2,2);plot(d{1});title('d')
%   suptitle(sprintf('Ex1, Identity, t=%g',t1))
%
% % Example 2: Linear Operator
%   nd=prod(size(m{im}));
%   forward.G=precal_cov_2d(nd,1,1,1,'1 Sph(10)');
%   tic;
%   d=sippi_forward(m,forward);
%   t2=toc;
%   figure(2);
%   subplot(1,2,1);plot(m{1});title('m')
%   subplot(1,2,2);plot(d{1});title('d')
%   suptitle(sprintf('Ex2, t=%g',t2))
%
%
%  % Example 3: Force forward.G to be sparse
%   forward.force_sparse=1;
%   tic;
%   d=sippi_forward(m,forward);
%   t3=toc;
%   figure(3);
%   subplot(1,2,1);plot(m{1});title('m')
%   subplot(1,2,2);plot(d{1});title('d')
%   suptitle(sprintf('Force sparse, t=%g',t3))
%
% See also: sippi_forward
function [d,forward,prior,data]=sippi_forward_linear(m,forward,prior,data,id,im)

if nargin<6,
    im=1;
end
if nargin<5,
    id=1;
end

%% setup forward structure
if ~isfield(forward,'is_identity')
    forward.is_identity=0;
end


if ~isfield(forward,'G')
    forward.G=speye(prod(size(m{im})));
    forward.is_identity=1;   
end


if ~isfield(forward,'force_sparse')
    forward.force_sparse=0;
end

if forward.force_sparse==1;
    if ~issparse(forward.G)
        forward.G=sparse(forward.G);
    end
end


%% solve forward problem
if forward.is_identity==1;
    d{id}=m{im}(:);
else
    d{id}=forward.G*m{im}(:);
end
