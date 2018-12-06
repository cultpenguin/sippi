% sippi_forward_linefit Line fit forward solver for SIPPI
%
%
% [d,forward,prior,data]=sippi_forward_fdem1d(m,forward,prior,data);
%
% forward.ds=0; % DOWNSAMPLING [1]:yes, [0]:no
% forward.S; % SYSTEM DESCRIPTION, see fdem1d
%
% forward.force_one_thread [0], 0:force one thread only ( nor parallellization)
% forward.nthreads : optionally maually choose the number of seperate
%                    Matlab workers for parfor
%

function [d,forward,prior,data]=sippi_forward_fdem1d(m,forward,prior,data);

if ~isfield(forward,'force_one_thread'); forward.force_one_thread=0;end
if ~isfield(forward,'ds'); forward.ds=0;end
if ~isfield(forward,'S');
    forward.S = readSystem('hemSystem.txt');
end
if ~isfield(forward,'ndata');
    forward.ndata=2*length(forward.S.freq);
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

%%
forward.model.CON=1./exp10(m{1}(:,forward.ix));
if (~isfield(forward.model,'X')|~isfield(forward.model,'Z'));
    [forward.model.X,forward.model.Z]=meshgrid(prior{1}.x(forward.ix),prior{1}.y);
end
if ~isfield(forward.model,'EL')
    forward.model.EL=repmat(prior{1}.y(:),[1 length(forward.x)]);
end

htx=-30;

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
        [d_est_par(:,i)] = calcFwd(forward.S,forward.model,htx,forward.ds,id_comp(i));
    end
    if length(id_comp)>0
        d_est(:,id_comp)=d_est_par;
    end
else
    for i=id_comp;%1:length(forward.x);
        [d_est(:,i)] = calcFwd(forward.S,forward.model,htx,forward.ds,i);
    end
end

d{1}=d_est(:);

