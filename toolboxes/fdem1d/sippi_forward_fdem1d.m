% sippi_forward_fdem1d: 1D fdem1d forward solver
%
%
% [d,forward,prior,data]=sippi_forward_fdem1d(m,forward,prior,data);
%
% forward.ds=0; % DOWNSAMPLING [1]:yes, [0]:no
% forward.S; % SYSTEM DESCRIPTION, see fdem1d
% forward.htx; % Height of TX below surface (negative above surface)
%
% forward.force_one_thread [0], 0:force one thread only ( nor parallellization)
% forward.nthreads : optionally maually choose the number of seperate
%                    Matlab workers   for parfor
%
%
% From; 
%  Elwaseif, M., Robinson, J., Day-Lewis, F. D., Ntarlagiannis, D., Slater, L. D., Lane, J. W., ... & Schultz, G. (2017). A matlab-based frequency-domain electromagnetic inversion code (FEMIC) with graphical user interface. Computers & Geosciences, 99, 61-71.
%  Minsley, B. J. (2011). A trans-dimensional Bayesian Markov chain Monte Carlo algorithm for model assessment using frequency-domain electromagnetic data. Geophysical Journal International, 187(1), 252-272.
%
% See also sippi_forward_fdem1d_example.m for an example.
%
%
%

function [d,forward,prior,data]=sippi_forward_fdem1d(m,forward,prior,data);

if ~isfield(forward,'force_one_thread'); forward.force_one_thread=0;end
if ~isfield(forward,'onedim'); forward.onedim=0;end
if ~isfield(forward,'ds'); forward.ds=0;end
if ~isfield(forward,'S');
    forward.S = readSystem('hemSystem.txt');
end
if ~isfield(forward,'ndata');
    forward.ndata=2*length(forward.S.freq);
end
if ~isfield(forward,'htx_as_data');
    forward.htx_as_data=0;
end
if (prior{1}.dim(1)>1)&&(prior{2}.dim(1)==1)
    % condunctivites  profiles in x direction
    forward.onedim=1;
    forward.x=[0];
else
    % condunctivites as columns in y direction
    forward.onedim=0;
    forward.x=prior{1}.x;    
end
    

% find index of dat locations
%% find location index of forward.x in prior model
if ~isfield(forward,'ix')    
    if (length(forward.x)==length(prior{1}.x))
        forward.ix=1:length(forward.x);
    else
        for ix=1:length(forward.x)
            iix=interp1(prior{1}.x,1:1:length(prior{1}.x),forward.x(ix),'nearest','extrap');
            iix=round(iix);
            forward.ix(ix)=iix;
        end
    end
end

forward.ix;

%%
if forward.onedim==1;
    forward.model.CON=1./exp10(m{1}(:));
    if (~isfield(forward.model,'X')|~isfield(forward.model,'Z'));
        [forward.model.X,forward.model.Z]=meshgrid(prior{1}.x,prior{1}.y);
        forward.model.Z=prior{1}.x(:);
        forward.model.X=0.*forward.model.Z;
    end
else
    forward.model.CON=1./exp10(m{1}(:,forward.ix));
    if (~isfield(forward.model,'X')|~isfield(forward.model,'Z'));
        [forward.model.X,forward.model.Z]=meshgrid(prior{1}.x(forward.ix),prior{1}.y);
    end
end
%if ~isfield(forward.model,'EL')
%    forward.model.EL=repmat(prior{1}.y(:),[1 length(forward.x)]);
%end

imm=length(prior);
for im=imm
    if strcmp(lower(prior{im}.name),'htx')
        forward.htx=-1*m{im};
    end
end

if ~isfield(forward,'htx')
    forward.htx=-30;
end
%disp(sprintf('htx=%5.2f',forward.htx));

%% PARALLEL
isOpen=0; % BY DEF, NO PARFOR/PARRALLELISATION
if forward.force_one_thread==0
    try % CHECK IF PARFOR IS AVAILABLE
        isOpen = isempty(gcp('nocreate'))==0;
        if ~isOpen
            sippi_verbose(sprintf('%s : SETTING PARRALLEL PARFOR LOOP',mfilename))
            myCluster = parcluster('local');
            if isfield(forward,'nthreads');
                myCluster.NumWorkers=forward.nthreads;
            end
            parpool(myCluster);
        end
        isOpen = ~isempty(gcp('nocreate'));
    end
end

d_est=zeros(forward.ndata,length(forward.x));

id_comp=1:length(forward.x);

if isfield(forward,'preferential_forward')
    if isfield(forward,'m_current');
        d_est=reshape(forward.d_current{1},12,length(forward.x));
        dm=mean(abs(m{1}-forward.m_current{1}));
        dm=dm(forward.ix); % consider only data for defines x locaitons in forward
        id_comp=find(dm>.00005);
    end   
end

%isOpen=0;
if (isOpen)
    parfor i=1:length(id_comp);%1:length(forward.x);
        [d_est_par(:,i)] = calcFwd(forward.S,forward.model,forward.htx,forward.ds,id_comp(i));
    end
    if length(id_comp)>0
        d_est(:,id_comp)=d_est_par;
    end
else
    for i=id_comp;%1:length(forward.x);
        [d_est(:,i)] = calcFwd(forward.S,forward.model,forward.htx,forward.ds,i);
    end
end

d{1}=d_est(:);
if forward.htx_as_data==1
    d{1}(end+1)=-1*forward.htx;
end



