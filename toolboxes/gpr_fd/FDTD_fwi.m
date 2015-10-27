
function [dt nt addpar]=FDTD_fwi(Eps,Sig,dx_fwd,t,ant_pos,addpar)

%================ FDTD simulation of electromagnetic waves ================
%
% Call: [dt nt]=FDTD_fwi(Eps,Sig,dx_fwd,t,ant_pos,addpar);
%
% General: The first five input parameters have to be defined and several additional 
% model parameters may be defined if need. The last additional 'addpar' input parameter 
% contains several parameters which may be defined if needed. If the additional 
% parameters are not defined default values are used.
%
% --- General input parameters (needed)---
% * Eps: 2D array containing the eps-model
% * Sig: 2D array containing the sig-model
% * dx_fwd: Cell size (m) for forward simulations
% * t: Total observation time (ns)
% * ant_pos: Nx4-array containing the (x,z)-tranitter and (x,z)-receiver positions.
%   N is the number of transmitter and receiver combination (Ntrn x Nrec)
%  
% --- General output parameters ----
% * nt: Number of samples in the output traces
% * dt: Sample interval of the output traces
%
% ------------------------------------------------------------------------
% The following parameters are only included if necessary, otherwise the
% default value is used. These additional parameters have to be defined
% as members in a structure: (addpar.(the exact parameter name as defined below))
% Ex: addpar.cores=2 will make the code run in multithreating mode on two cores. 
%
% --- Executables ---
% * forwardexeN: Contains strings with the names of the executables used 
%   for the forward simulation. N is an integer. N equals the number of 
%   cores applied during the forward simulation. If one core is applied 
%   only the file forwardexe1 is need. If multiple cores are used for 
%   the simulation more forwardexe files have to be defined.  
%
% --- Executable --- [NEW]
% addpar.executable complete path to executable
%      default directory is the directory in which FDTD_fwi resides
%      default executable name is 'FDTD_forward' appended with the
%      architecture as output from "computer('arch')". such as for example
%      "FDTD_forward_glnxa64" on 64 bit linux and 
%      "FDTD_forward_maci64" on 64 bit OSX
%      in windows the executable is simply
%      "FDTD_forward.exe"
%
% --- Output working directory ---
%  addpar.work_dir
%
% --- Boundaries ---
% * bx: The number of boundary cells (default=40)
% * sm: GPML-Tuning; Constant characterising the coordinate stretching in
%   the boundaries (-1:auto) (default=auto)
% * sigm: GPML-Tuning; Constant describing the maximum of the conductivity
%   profile in the boundaries (-1:auto) (default=auto)
% * n: GPML-Tuning; Exponent, describing how fast the stretching and
%   conductivity happens in the boundaries (default=2)
%
% --- Stabilisation ---
% * sc: Courant stability scaling (default=0.9)
% * Epsmin: Smallest eps value used to calculate the time increment used in the 
%   simulation (in units of relative dielectric permittivity)
%   (-1:auto) (default=1)
%
% --- Source wavelet ---
% * srcWType: Source-wavelet-type (1:Gauss, 2:Diff. Gauss, 3:Ricker,
%   4:External) (default=Gauss)
% * srcImpl: Source implementation; (0:hard-source) or (1:soft-source, that
%   accounts for the field at source-position) (default=soft-source)
% * srcOri: Source orientation (0:Ez) or (1:Er/Ex) (default=Ez)
% * Tg: Total Pulse Length (in s) if source-type=1 or 2 (default=1.4*10^-8 s);
%   Center frequency (in Hz) if source-type=3 (default=14*10^6 Hz)
% * Tp: Pulse Width (in s) only if source-type=1 or 2; Width of the pulse at
%   amplitude=exp(-0.5)*ampl_max (default=1.6*10^-9 s)
%
% --- Magnetic properties ---
% * Mu: Magnetic Permeability, which is assumed constant across the model
%   (in units of magnetic permeability of free space) (default=1)
% * rho: Equivalent magnetic loss or "magnetic resistivity" (default=0)
%
% --- Coordinate system ---
% * coord: Coordinate system in which the simulation is performed 
%   (0:Cartesian, 1:Cylindrical) (default=Cartesian)
% * dx_inv: The cell size for the inversion parameters. Is used to determine
%   where to position the receivers for the forward-backward simulation.
%
% --- Multithreading ----
% * cores: Number of cores applied for the forward simulation (default=1).
%   The number of cores will automatically be reduced if the number of chosen 
%   cores exceeds the number of transmitters. 
%
% -- Output --
% * snap: Snapshot "frequency"; Number of time-steps between snapshots 
%   saved to disk (default=10000 (~no snapshot output)) 
% * frac: Factor by which snapshots are resampled (to save disk space) (default=1)
% * sim_mode: Types of output used for pure forward (=1), forward
%   simulation followed by a backword simulation used for inversion purpose
%   (=2) and pure forward simulation for inversion purpose (=3) (e.g. forward
%   simualtion in a perturbed model) (default=1).
% * output: Determines which field direction to be outputted; (1)Ez, 
%   (2)Ex/Er, (3)Both Ez and Ex/Er, (4) The same as (1) but the numbering is 
%   given by the vector "no". This option is only active when sim_mode=1 is 
%   applied (default=1). 
% * no: is a vector containing numbers used for the numbering of the output
%   files. length(no) has to equal the number of transmitters.
% * debug: 1 = write parameter information to screen. 0 = nothing is
%   written to the screen.
%
% --- Various ---
% * start: Flag controlling if the simulation is executed or only the
%   parameter files are written (1:"on", 0:"off") (default="off") 
% * plot: Flag controlling if a plot of the setup it imaged (1:"plot", 0:"no plot") 
%   (default="no plot")
%
%
% The full-waveform FDTD forward modelling program executed by this function is 
% written by Jacques Ernst (ETH Zurich, 2006) as part of his Ph.D. project.
%
% Knud Cordua (2008-2015)                                            v. 1.2
% Thomas Mejer Hansen (2015)
%==========================================================================

% Debuging:
try 
    isempty(addpar.debug);
catch
    addpar.debug=0;
end

% executables and working directory
if ~isfield(addpar,'executable');
  root_dir=fileparts(which('FDTD_fwi.m'));
  if isunix
    f_executable=[root_dir,filesep,'FDTD_forward_',computer('arch')];
    if ~exist(f_executable,'file')
      f_executable=[root_dir,filesep,'FDTD_forward'];
    end
    addpar.executable=f_executable;
  else
    addpar.executable=[root_dir,filesep,'FDTD_forward.exe'];
  end
  if addpar.debug>-1
    disp(sprintf('Using binary: %s ',addpar.executable));
  end
end


if ~isfield(addpar,'workdir');
  addpar.workdir=pwd;
  if addpar.debug>-1
    disp(sprintf('Using working directory: %s ',addpar.workdir));
  end
end

% Forward and inverse cell size in m and simulation time in ns:
addpar.dx=dx_fwd;
addpar.t=t;

% Determine the size of the models in cells:
[addpar.Epsz,addpar.Epsx]=size(Eps);
[addpar.Sigz,addpar.Sigx]=size(Sig);

% Default values in case of error meassages:
dt=NaN;
nt=NaN;

% In case no additional parameters are chosen:
if nargin<6
    disp('-----------------------------------------------')
    disp('Only default values are applied')
end

% Default values in case they are not assinged later on (e.g. addpar.start=0)
dt=NaN;
nt=NaN;

% Compare the Sig and Eps model:
if addpar.Epsz~=addpar.Sigz
    disp('Warning: The depth of the Eps model is different from the depth of the Sig model')
    disp('-----------------------------------------------')
    return
else
    addpar.nzce=addpar.Epsz;
end

% Number of elements in x-direction:
if addpar.Epsx~=addpar.Sigx
    disp('Warning: The width of the Eps model is different from the width of the Sig model')
    disp('-----------------------------------------------')
    return
else
    addpar.nxce=addpar.Epsx;
end

% Write the chosen input parameters to the screen
if addpar.debug==1
    [addpar error]=write_parameters_to_screen(ant_pos,Sig,Eps,addpar);
    if error==1
    disp('Warning: sumthing went wrong in the function "write_parameters_to_screen"');
    return
    end
else
    [addpar error]=setup_input_parameters(ant_pos,Sig,Eps,addpar);
    if error==1
    disp('Warning: sumthing went wrong in the function "setup_input_parameters"');
    return
    end
end

% Check if any transmitter or receiver-position is outside the border of the model
I=find(ant_pos(:,1)>addpar.Epsx*addpar.dx);
if ~isempty(I)
    disp('Warning: A transmitter-position has exceeded the x-axis')
    disp('-----------------------------------------------')
    return
end
I=find(ant_pos(:,2)>addpar.Epsz*addpar.dx);
if ~isempty(I)
    disp('Warning: A transmitter-position has exceeded the z-axis')
    disp('-----------------------------------------------')
    return
end
I=find(ant_pos(:,3)>addpar.Epsx*addpar.dx);
if ~isempty(I)
    disp('Warning: A receiver-position has exceeded the x-axis')
    disp('-----------------------------------------------')
    return
end
I=find(ant_pos(:,4)>addpar.Epsz*addpar.dx);
if ~isempty(I)
    disp('Warning: A receiver-position has exceeded the z-axis')
    disp('-----------------------------------------------')
    return
end
I=find(ant_pos(:,2)<addpar.dx/2);
if ~isempty(I)
    disp('Warning: A transmitter-position is too close to the top')
    disp('-----------------------------------------------')
    return
end
I=find(ant_pos(:,4)<addpar.dx/2);
if ~isempty(I)
    disp('Warning: A receiver-position is too close to the top')
    disp('-----------------------------------------------')
    return
end


% Distribute the parameter files to the number of cores applied
[Ncores_applied error]=distribute_to_cores(ant_pos,Eps,Sig,addpar);
if error==1
    disp('Warning: Something went wrong in the subfunction "distribute_to_cores"')
    return
end
if addpar.debug==1
disp('-----------------------------------------------')
end

% Starting the forward simulation
if addpar.start==1
    [dt nt error]=fwi_execute(Ncores_applied,ant_pos,addpar.sim_mode,addpar);
    if error==1
        disp('Something went wrong in the subfunction "fwi_execute"')
        return
    end
end

if addpar.start==1 && addpar.debug==1
% Display some of the results from the forward simulation
disp(sprintf('%s %4.0f','Number of time-steps:',nt))
disp(sprintf('%s %e %s','Time increment:',dt*10^9,'ns'))
disp('-----------------------------------------------')
end