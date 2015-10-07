
function [Ncores_applied error]=distribute_to_cores(ant_pos,Eps,Sig,addpar)

error=0;
Ncores=addpar.cores;

% Determine the number of transmitter positions
try
    Ntrn=1;
    pos_new_trn(Ntrn)=1;
    for i=2:length(ant_pos(:,1))
        if ant_pos(i,1)~=ant_pos(i-1,1) | ant_pos(i,2)~=ant_pos(i-1,2)
            Ntrn=Ntrn+1;
            pos_new_trn(Ntrn)=i;
        end
    end
catch
    error=1;
    return
end
if addpar.debug==1
disp(sprintf('Number of transmitters: %i',Ntrn))
end


% Determine the number of transmitters to be simulated on each core
N_trn_pr_core=floor(Ntrn/Ncores);
Core_nfo(1:Ncores)=N_trn_pr_core;
if N_trn_pr_core*Ncores~=Ntrn
    Remain_trn=Ntrn-N_trn_pr_core*Ncores;
    for i=1:Remain_trn
        Core_nfo(i)=Core_nfo(i)+1;
    end
end

% Determine if the #trn<#cores
if length(Core_nfo)>sum(Core_nfo)
    for i=1:length(Core_nfo)
        if Core_nfo(i)~=0
            Core_nfo_tmp(i)=Core_nfo(i);
        end
    end
    Core_nfo=Core_nfo_tmp;
    Ncores=length(Core_nfo);
    if addpar.debug==1
    disp(sprintf('Core reduction: Only %i cores are applied (Reason: #trn<#cores)',Ncores))
    end
end

Ncores_applied=0;
Count_trn_old=1;
if ~exist(addpar.workdir,'dir');
  try
    mkdir(addpar.workdir);
  catch
    disp(sprintf('%s: could name make directory %s',mfilename,addpar.workdir))
  end
end

%%
for i=1:Ncores
    %
    WorkDir=[addpar.workdir,filesep,sprintf('run%03d',i)];
    if ~exist(WorkDir,'dir');
      mkdir(WorkDir);
    end

    % Find the path of the forward executable:
    %fpath=fileparts(which(eval(sprintf('%s%i','addpar.forwardexe',i))));
    fpath=WorkDir;
    
    % Delete old files before new files are generated
    delete(sprintf('%s',fpath,filesep,'*.fwd'))
    delete(sprintf('%s',fpath,filesep,'*.eps'))
    delete(sprintf('%s',fpath,filesep,'*.sig'))
    delete(sprintf('%s',fpath,filesep,'*.inv'))
    
    if i<Ncores && sum(Core_nfo)>1
        Count_trn=sum(Core_nfo(1:i))+1;
        tmp_pos=ant_pos(pos_new_trn(Count_trn_old):pos_new_trn(Count_trn)-1,:);
        Count_trn_old=Count_trn;
        Ncores_applied=Ncores_applied+1;
    else
        tmp_pos=ant_pos(pos_new_trn(Count_trn_old):end,:);
        Ncores_applied=Ncores_applied+1;
    end
    
    
    %========= Generating the input parameters for simulation ================%
    %======================= source.E=======================%
    if isfield(addpar,'wavelet');
      save_wavelet(addpar.wavelet.data,addpar.wavelet.dt,[fpath,filesep,'source.E']);
    end
    
    %======================= Create radar.FWD File =======================%
    f_radar_fwd=sprintf('%s%s',fpath,filesep,'radar.fwd');
    Fid=fopen(f_radar_fwd, 'wt');
    fprintf ( Fid,'----------------------- Radar-Forward-Config ------------------------\n');
    fprintf ( Fid,'File-Version\t\t\t\t=%i\n',                      7);
    fprintf ( Fid,'INV-FWD-Scale Factor\t\t\t=%i\n',                1);
    fprintf ( Fid,'#Boundary-Cells\t\t\t\t=%i\n',                   addpar.bx);
    %use for example (lamda_min/(1+sm)>2* to 3*dx) -> 0=simple wave prop.
    fprintf ( Fid,'GPML-Tuning: sm (-1:auto)\t\t=%f\n',             addpar.sm);
    %use for example ([(eps*c/(n_bou*dx/InvFwd_F))/(1+sm*(1/3+2/pi^2))]*ln(R_th))
    fprintf ( Fid,'GPML-Tuning: sigm (-1:auto)\t\t=%f\n',           addpar.sigm);
    fprintf ( Fid,'GPML-Tuning: expo\t\t\t=%f\n',                   addpar.n);
    fprintf ( Fid,'Total Observation Time\t\t\t=%3.2e\n',           addpar.t);
    fprintf ( Fid,'C-Stability-Scaling\t\t\t=%f\n',                 addpar.sc);
    fprintf ( Fid,'Min. Eps (Stab.; Default:1;auto:-1)\t=%4.2f\n',  addpar.Epsmin);
    %Source Parameters
    fprintf ( Fid,'Source-Type (1:G, 2:DG, 3:R, 4:Ex)\t=%i\n',      addpar.srcWType);
    fprintf ( Fid,'Inc/Exc Field at SRC\t\t\t=%i\n',                addpar.srcImpl);
    fprintf ( Fid,'SRC orientation (0:Ez, 1:Er/Ex)\t\t=%i\n',       addpar.srcOri);
    switch addpar.srcWType
        case { 1,2 }
            fprintf ( Fid,'Total Pulse Length\t\t\t=%3.2e\n',       addpar.Tg); %*1e-9);
            fprintf ( Fid,'Pulse Width\t\t\t\t=%3.2e\n',            addpar.Tp); %*1e-9);
        case { 3 }
            fprintf ( Fid,'Center Frequency\t\t\t=%3.2e\n',         addpar.Tg); %*1e6);
            fprintf ( Fid,'Pulse Width\t\t\t\t=%i\n',               -1);
        case { 4 }
            fprintf ( Fid,'Total Pulse Length\t\t\t=%i\n',          -1);
            fprintf ( Fid,'Pulse Width\t\t\t\t=%i\n',               -1);
    end;
    fprintf ( Fid,'Magnetic Permeability\t\t\t=%4.2f\t\n',          addpar.mu); %in MU0!!!!
    fprintf ( Fid,'Equivalent Magnetic Loss\t\t=%4.2f\n',           addpar.rho);
    fclose (Fid);

    %======================= Create radar.INV File =======================%
    f_radar_inv=sprintf('%s%s%s',fpath,filesep,'radar.inv');
    Fid=fopen(f_radar_inv, 'wt');
    fprintf ( Fid,'---------------------- Radar-Inversion-Config -----------------------\n');
    fprintf ( Fid,'File-Version\t\t\t\t=%i\n',                      7);
    fprintf ( Fid,'DATA-File Path\t\t\t\t=%s\\\n',                  '\');
    fprintf ( Fid,'Output Path(Disabled)\t\t\t=%s\\\n',             '\');
    fprintf ( Fid,'Small Perturbation (StepLCalc)\t\t=%3.2e\n',     1.00e-5);
    fprintf ( Fid,'CAUTION: Gradient Scaling\t\t=%3.2e\n',          1);
    fprintf ( Fid,'CAUTION: Disable Trc Norm. (0:Off,1:On)\t=%i\n', 0);
    fprintf ( Fid,'Bleed Zone (0:Off,1:On)\t\t\t=%i\n',             1);
    fprintf ( Fid,'Gradient Smoothing (0:Off,1:On)\t\t=%i\n',       1);
    fprintf ( Fid,'-> GS Filter (1:sG,2:fG)\t\t=%i\n',              0);
    fprintf ( Fid,'Trn/Rec Muting (0:Off,1:On)\t\t=%i\n',           0);
    fprintf ( Fid,'-> T/R Muting: R(tot)\t\t\t=%i\n',               0);
    fprintf ( Fid,'-> T/R Muting: R(zero)\t\t\t=%i\n',              0);
    fprintf ( Fid,'Trace Marker File (0:off,1:on)\t\t=%i\n',        0);
    fprintf ( Fid,'Trc. Mute: Energy Threshold\t\t=%3.1f\n',        0);
    fprintf ( Fid,'Eps(1) or Sig(2) Inversion\t\t=%i\n',            1);
    fprintf ( Fid,'SRC Inversion (0:Off,1:On)\t\t=%i\n',            0);
    fprintf ( Fid,'Start Iteration Nr.\t\t\t=%i\n',                 0);
    fprintf ( Fid,'Total #Iterations\t\t\t=%i\n',                   0);
    fprintf ( Fid,'Cell-Size\t\t\t\t=%5.3f\n',                      addpar.dx);
    fclose ( Fid);

    %======================= Create model.eps File =======================%
    save2Dbin(sprintf('%s%s%s',fpath,filesep,'model.eps'), addpar.nxce, addpar.nzce, Eps);

    %====================== Create model.sig File ========================%
    save2Dbin(sprintf('%s%s%s',fpath,filesep,'model.sig'), addpar.nxce, addpar.nzce, Sig);
end

% ============= Function used to write the binary model files ============%

function save2Dbin(fname,nx,nz,data)
Fid=fopen(sprintf('%s',fname),'wb');
fwrite(Fid,nx,'int');
fwrite(Fid,nz,'int');
fwrite(Fid,data,'double');
fclose(Fid);

