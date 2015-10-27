% sippi_forward_gpr_fd: full waveform gpr forward
%
% Call :
%  [d,forward,prior,data]=sippi_forward_gpr_fd(m,forward,prior,data,id,im)
%
%  the prior must be such that m{1} relfect the eps field, and (optionally)
%  m{2} reflect the sig field (if not set it is trated as constance).
% 
% % Mandatory 
% forward.ant_pos: antenna postions for receivers and recorders
%                  [Sx1, Sy1, Rx1, Ry1    
%                   Sx2, Sy2, Rx2, Ry2
%                   ..];
%
% % Optional settings
% forward.T=100*10^-9; % Frequecy of sourve wavelet
% forward.sig=3; % If not set it is set as a constant field
% forward.output_type='shot'; % each shot gather is output as individual data structures
% forward.output_type='trace'; % each trace is output as individual data structures
% forward.output_it=10; % only output every 'output_it' time sample
%
%
% forward.dx_forward % the spatial sampling used for FDTD modeling
%                      If this is different than the spatial sampling
%                      density for the prior, then prior realizations are
%                      rescales before forward modeling
%
%
% % Settings for forward modeler
% All seetings for the forward modelling code can be set in the 
% forward.addpar structure. E.g to set the number of threads used:
% forward.addpar.cores=4;%
%
% See also FDTD_fwi
%
%
function [d,forward,prior,data]=sippi_forward_gpr_fd(m,forward,prior,data,id,im)

if nargin<6, im=1;end
if nargin<5, id=1;end
if nargin<4, data{1}.null='';end

eps=m{1};
if length(m)>1
  sig=m{2};
else
  if ~isfield(forward,'sig');
    forward.sig=3;
  end
  sig=ones(size(m{1})).*forward.sig;
end

% output only every 'output_it' sample
if ~isfield(forward,'output_it');
  forward.it=1;
end

% return in data in structures of SHOT of TRACE data
if ~isfield(forward,'output_type');
  %forward.output_type='shot';
  forward.output_type='trace';
end
  
if ~isfield(forward,'dx_fwd');
    dx_in=prior{1}.x(2)-prior{1}.x(1);
    forward.dx_fwd=dx_in;
end


if ~isfield(forward,'t');
  forward.t=1e-7;
end

if ~isfield(forward,'EPS0');
  forward.EPS0=8.85418781762039080*1e-12;
end
if ~isfield(forward,'SIG0');
  forward.SIG0=10^-3;
end
if ~isfield(forward,'MU0');
  forward.MU0=1.25663706143591730*1e-6; %Magnetic permeability of free space, Vs/Am (4*pi*10^-7 N/A^2);
end
if ~isfield(forward,'c0');
  forward.c0=2.99792458*10^8;
end

if ~isfield(forward,'addpar'); forward.addpar.null=[];end
if ~isfield(forward.addpar,'cores');
    % try to locate the number of logical cores
    try
      ncores=java.lang.Runtime.getRuntime().availableProcessors;
    catch
      ncores=1;
    end
    forward.addpar.cores=ncores;
end
if ~isfield(forward.addpar,'status');forward.addpar.status=1;end
if ~isfield(forward.addpar,'plot');forward.addpar.plot=0;end
if ~isfield(forward.addpar,'srcWType');forward.addpar.srcWType=3;end
if ~isfield(forward.addpar,'Tg');
  if isfield(forward,'Tg')
    forward.addpar.Tg=forward.Tg;  
  else
    forward.addpar.Tg=100*10^6;
  end
end
if ~isfield(forward.addpar,'Epsmin');forward.addpar.Epsmin=1;end
if ~isfield(forward.addpar,'start');forward.addpar.start=1;end
if ~isfield(forward.addpar,'debug');forward.addpar.debug=-1;end

%% run forward
dx_in=prior{1}.x(2)-prior{1}.x(1);
if dx_in~=forward.dx_fwd;
    eps_fd=resize(eps,dx_in/forward.dx_fwd)*forward.EPS0;
    sig_fd=resize(sig,dx_in/forward.dx_fwd)*forward.SIG0;
else
    eps_fd=eps.*forward.EPS0;
    sig_fd=sig.*forward.SIG0;
end

[forward.dt forward.nt forward.addpar]=FDTD_fwi(eps_fd,sig_fd,forward.dx_fwd,forward.t,forward.ant_pos,forward.addpar);

forward.time=[0:(forward.nt-1)].*forward.dt;

%% get data
% get output
tpos=unique(forward.ant_pos(:,1:2),'rows');
forward.Nsources=size(tpos,1);

path_sim=[pwd,filesep,'observed_Ez_trn'];


j=0; % counter for trace id
for id=1:forward.Nsources
  filename_sim=sprintf('%s%d%s',path_sim,id,'.mat');
  try
    load(filename_sim);
  catch
    keyboard
  end
  d_sim=Trace_Ez';
  d_sim=d_sim(forward.output_it:forward.output_it:end,:);
  if strcmp(forward.output_type,'shot')
    d{id}=d_sim(:);
  elseif strcmp(forward.output_type,'trace')
    for it=1:size(d_sim,2)
      j=j+1;
      d{j}=d_sim(:,it);
    end
  end
  
end

