calc\_Cd {#calc_Cd}
--------

      Calc_cd Setup a covariance model to account for borehole imperfections
     
      Call: Cd=calc_Cd(ant_pos,var_uncor,var_cor1,var_cor2,L)
      This function sets up a data covariance matrix that accounts for static
      (i.e. correlated) data errors.
      
      Inputs:
      * ant_pos: A N x 4 array that contains N combinations of transmitter/source 
      and receiver positions. The first two columns are the x- and y-coordinates
      of the transmitter/source position. The last two columns are the x- and 
      y-coordiantes of the receiver position.
      * var_uncor: The variance of the uncorrelated data errors.
      * var_cor1: The variance of the correlated data errors
      related to the transmitter/source positions.
      * var_cor2: The variance of the correlated data errors
      related to the receiver positions.
      * L: The correlation length for the correlation between the individual
      transmitter/source or receiver positions using an exponential covariance 
      function. For typical static errors the correlation length is set to a 
      small number (e.g. 10^-6).
      
      For more details and practical examples see:
      Cordua et al., 2008 in Vadose zone journal.
      Cordua et al., 2009 in Journal of applied geophysics.
     
      Knud S. Cordua (2012)

eikonal
-------

      eikonal Traveltime computation by solving the eikonal equation
      
      tmap=eikonal(x,y,z,V,Sources,type);
     
       x,y,z : arrays defining the x, y, and z axis
       V: velocity field, with size (length(y),length(x),length(z));
       Sources [ndata,ndim] : Source positions
       type (optional): type of eikonal solver: [1]:Fast Marching(default), [2]:FD
     
       tmap [size(V)]: travel times computed everywhere in the velocity grid
     
     %Example (2D):
        x=[1:1:100];
        y=1:1:100;
        z=1;
        V=ones(100,100);V(:,1:50)=2;
        Sources = [10 50;75 50];
        t=eikonal(x,y,z,V,Sources);
        subplot(1,2,1);imagesc(x,y,t(:,:,1,1));axis image;colorbar
        subplot(1,2,2);imagesc(x,y,t(:,:,1,2));axis image;colorbar
     
      See also eikonal_traveltime
     

eikonal\_raylength
------------------

      eikonal_raylength : Computes the raylength from S to R using the eikonal equaiton
     
      Call:
        raylength=eikonal_raylength(x,y,v,S,R,tS,doPlot)
     

eikonal\_traveltime
-------------------

      eikonal_traveltime Computes traveltime between sources and receivers by solving the eikonal equation
     
      t=eikonal_traveltime(x,y,z,V,Sources,Receivers,iuse,type);
     
       x,y,z : arrays defining the x, y, and z axis
       V: velocity field, with size (length(y),length(x),length(z));
       Sources [ndata,ndim] : Source positions
       Receivers [ndata,ndim] : Receiver positions
       iuse (optional): optionally only use subset of data. eg.g i_use=[1 2 4];
       type (optional): type of eikonal solver: [1]:Fast Marching(default), [2]:FD
     
       tmap [size(V)]: travel times computed everywhere in the velocity grid
     
     %Example (2%
     
      Example 2d traveltime compuation
     
      Example (2D):
        x=[1:1:100];
        y=1:1:100;
        z=1;
        V=ones(100,100);V(:,1:50)=2;
        S=[50 50 1;50 50 1];
        R=[90 90 1; 90 80 1];
        t=eikonal_traveltime(x,y,z,V,S,R)
     
      Example (3D):
        nx=50;ny=50;nz=50;
        x=1:1:nx;
        y=1:1:ny;
        z=1:1:nz;
        V=ones(ny,nx,nz);V(:,1:50,:)=2;
        S=[10 10 1;10 10 1;10 9 1];
        R=[40 40 40; 40 39 40; 40 40 40];
        t=eikonal_traveltime(x,y,z,V,S,R)
     
     
      See also eikonal

kernel\_buursink\_2d
--------------------

      kernel_buursink_2k Computes 2D Sensitivity kernel based on 1st order EM scattering theory
      
      See 
        Buursink et al. 2008. Crosshole radar velocity tomography 
                              with finite-frequency Fresnel. Geophys J. Int.
                              (172) 117;
     
       CALL : 
          % specify a source trace (dt, wf_trace):
          [kernel,L,L1_all,L2_all]=kernel_buursink_2d(model,x,z,S,R,dt,wf_trace);
          % Use a ricker wavelet with center frequency 'f0'
          [kernel,L,L1_all,L2_all]=kernel_buursink_2d(model,x,z,S,R,f0));
     
     
      Knud Cordua, 2009, 
      Thomas Mejer Hansen (small edits, 2009)
     

kernel\_finite\_2d
------------------

      kernel_finite_2d 2D sensitivity kernels
     
       Call:
         [Knorm,K,dt,options]=kernel_finite_2d(v_ref,x,y,S,R,freq,options);

kernel\_fresnel\_2d
-------------------

      kernel_fresnel_2d Sensitivity kernel for amplitude and first arrival
     
      Call:
        [kernel_t,kernel_a,P_omega,omega]=kernel_fresnel_2d(v,x,y,S,R,omega,P_omega);
     
     
      Based on Liu, Dong, Wang, Zhu and Ma, 2009, Sensitivity kernels for
      seismic Fresenl volume Tomography, Geophysics, 75(5), U35-U46
     
      See also kernel_fresnel_monochrome_2d
     
      Run with no argument for an example.
     
     

kernel\_fresnel\_monochrome\_2d
-------------------------------

      kernel_fresnel_monochrome_2d 2D monchrome kernel for amplitude and first arrival
     
      Call:
        [kernel_t,kernel_a]=kernel_fresnel_monochrome_2d(v,x,y,S,R,omega);
      or
        [kernel_t,kernel_a]=kernel_fresnel_monochrome_2d(v,x,y,S,R,omega,L,L1,L2);
     
      Based on Liu, Dong, Wang, Zhu and Ma, 2009, Sensitivity kernels for
      seismic Fresenl volume Tomography, Geophysics, 75(5), U35-U46
     
      See also, kernel_fresnel_2d
     

kernel\_multiple
----------------

      kernel_multiple Computes the sensitivity kernel for a wave traveling
      from S to R.
     
      CALL : 
         [K,RAY,Gk,Gray,timeS,timeR,raypath]=kernel_multiple(Vel,x,y,z,S,R,T,alpha,Knorm);
     
      IN : 
         Vel [ny,nx] : Velocity field
         x [1:nx] :
         y [1:ny] :
         z [1:nz] :
         S [1,3] : Location of Source
         R [1,3] : Location of Receiver
         T : Donminant period
         alpha: controls exponential decay away ray path
         Knorm [1] : normaliztion of K [0]:none, K:[1]:vertical
     
      OUT :
         K : Sensitivity kernel
         R : Ray sensitivity kernel (High Frequency approx)
         timeS : travel computed form Source
         timeR : travel computed form Receiver
         raypath [nraydata,ndim] : the center of the raypath 
     
      The sensitivity is the length travelled in each cell.
      
     
      See also : fast_fd_2d
     
      TMH/2006
     

kernel\_slowness\_to\_velocity
------------------------------

      kernel_slowness_to_velocity Converts from slowness to velocity parameterizations
     
      G : kernel [1,nkernels]
      V : Velocity field (
     
     
      CALL:
        G_vel=kernel_slowness_to_velocity(G,V);
      or 
        [G_vel,v_obs]=kernel_slowness_to_velocity(G,V,t);
      or
        [G_vel,v_obs,Cd_v]=kernel_slowness_to_velocity(G,V,t,Cd);
     
      

mspectrum
---------

      mspectrum : Amplitude and Power spectrum
      Call  :
          function [A,P,smoothP,kx]=mspectrum(x,dx)
      
      1D (A)mplitude and (P)owerspectrum of x-series with spacing dx
     

munk\_fresnel\_2d
-----------------

      2D frechet kernel, First Fresnel Zone 
     
      See Jensen, Jacobsen, Christensen-Dalsgaard (2000) Solar Physics 192.
     
      Call :
      S=munk_fresnel_2d(T,dt,alpha,As,Ar,K);
     
      T : dominant period
      dt : 
      alpha : degree of cancellation 
      As : Amplitude fo the wavefield propagating from the source
      Ar : Amplitude fo the wavefield propagating from the receiver
      K : normalization factor

munk\_fresnel\_3d
-----------------

      3D frechet kernel, First Fresnel Zone 
     
      See Jensen, Jacobsen, Christensen-Dalsgaard (2000) Solar Physics 192.
     
      Call :

pick\_first\_arrival
--------------------

      pick_first_arrival : pick first arrival travel time data using simple
                           correlation
     
      Call 
        [tt_pick]=pick_first_arrival(wf_data,ref_trace,ref_t0,doPlot,wf_time);
     

plot\_traveltime\_sr
--------------------

      plot_traveltime_sr
     
      Call 
         plot_traveltime(S,R)
         S: [n,2] : source locattion
         R: [n,2] : reveiver locattion
         or (3d)
         S: [n,3] : source locattion
         R: [n,3] : reveiver locattion
     
     
      EX:
       % 2D
       D=load('AM13_data.mat');
       plot_traveltime_sr(D.S,D.R);
       or
       ant_pos=[D.S,D.R];
       plot_traveltime_sr(ant_pos);
     
       % 3D
       D=load('AM1234_data.mat');
       plot_traveltime_sr(D.S,D.R);
     
     

sippi\_forward\_traveltime
--------------------------

      sippi_forward_traveltime Traveltime computation in SIPPI
     
      Call :
        [d,forward,prior,data]=sippi_forward_traveltime(m,forward,prior,data)
     
        forward.type determines the method used to compute travel times
        forward.type='ray';
        forward.type='fat';
        forward.type='eikonal';
        forward.type='born';
        forward.type='fd';
     
        forward.sources [ndata,ndim]: Source locations
        forward.receivers [ndata,ndim]: Receiver locations
     
        % the following options does not apply to 'eikonal' type modeling
        forward.linear : [0] a linear kernel is computed, based on the current velocity model
                         [1] a linear kenrel is computed only once, based on
                         the velocity field defined in forward.linear_m;
     
        forward.linear_m: the reference velocity field, for a linear forward
                          operator (forward.G) will be computed.
                          Can be eithe a scalar (constant velocity field) or
                          the same size as the the velcity model 'm'.
     
        forward.normalize_vertical [1]: Normalize sensitivitykernel by
                                        in vertical slices
                                   [0]: No normalization
     
        forward.alpha [1]: alpha value for munk_fresnel_2d, munk_fresnel_3f

sippi\_forward\_traveltime\_example
-----------------------------------

      sippi_forward_traveltime_example: Different examples of travel time
           computation
     
      See also: sippi_forward_traveltime
     
     % load some data

tomography\_kernel
------------------

      tomography_kernel Computes the sensitivity kernel for a wave traveling from S to R.
     
      CALL :
         [K,RAY,Gk,Gray,timeS,timeR,raypath]=tomography_kernel(Vel,x,y,z,S,R,T,alpha,Knorm);
     
      IN :
         Vel [ny,nx] : Velocity field
         x [1:nx] :
         y [1:ny] :
         z [1:nz] :
         S [1,3] : Location of Source
         R [1,3] : Location of Receiver
         T : Donminant period
         alpha: controls exponential decay away ray path
         Knorm [1] : normaliztion of K [0]:none, K:[1]:vertical
     
      OUT :
         K : Sensitivity kernel
         R : Ray sensitivity kernel (High Frequency approx)
         timeS : travel computed form Source
         timeR : travel computed form Receiver
         raypath [nraydata,ndim] : the center of the raypath
     
      The sensitivity is the length travelled in each cell.
     
     
     
