% sippi_forward_dcfd2_5D
%
% Call :
%   [d,forward,prior,data]=sippi_forward_dcfd2_5D(m,forward,prior,data)
%
%
function [d,forward,prior,data]=sippi_forward_dcfd2_5D(m,forward,prior,data,id,im)

if nargin<2
    forward.null=[];
end

if nargin<4;    data{1}.null='';end
if nargin<5;    id=1;end
if nargin<6;
    if isfield(forward,'im')
        im=forward.im;
    else
        im=1;
    end
end


% Test if Para has been defined
if ~isfield(forward,'Para');
    
    if ~isfield(forward,'srcloc')
        tmp=load('srcloc_hvede_unique_e257.mat');
        forward.srcloc=tmp.srcloc_hvede_unique;
    end
    
    if ~isfield(forward,'dx')
        tmp=load('dxdz_hvede_ai2.mat');
        forward.dx = tmp.dx_real_1;
    end
    if ~isfield(forward,'dz')
        tmp=load('dxdz_hvede_ai2.mat');
        forward.dz = tmp.dz_real_1;
    end
    
    if ~isfield(forward,'BC_cor')
        forward.BC_cor = [];
    end
    
    
    if ~isfield(forward,'num')
        forward.num = 4;
    end
    
    if ~isfield(forward,'recloc')
        tmp=load('recloc_hvede_e257.mat');
        forward.recloc = tmp.recloc_hvede;
    end
    
    if ~isfield(forward,'srcnum')
        tmp=load('srcnum_hvede_e257.mat');
        forward.srcnum = tmp.srcnum_hvede;
    end
    
    fprintf('\b%s: dcfd2_5F: Setting up Parameterization\n',mfilename)
    forward.Para = get_2_5Dpara(forward.srcloc,forward.dx,forward.dz,forward.BC_cor,forward.num,forward.recloc,forward.srcnum);
end

%load srcloc_hvede_unique_e257.mat % shifted laterally -3
%load dxdz_hvede_ai2.mat % final model from inversion, modif on upper cells
%load recloc_hvede_e257.mat % shifted laterally -3
%load srcnum_hvede_e257.mat
%Para=make_Para(srcloc_hvede_unique,dx_real_1,dz_real_1,recloc_hvede,srcnum_hvede);
%load s_hvede_ai.mat
%d{id} = dcfw2_5D(s_hvede_3,forward.Para);

d{id} = dcfw2_5D(m{1},forward.Para);


