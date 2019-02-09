sippi\_adjust\_step\_size
-------------------------

      sippi_adjust_step_size Adjust step length length for Metropolis sampler in SIPPI
       
      Call : 
        step=sippi_adjust_step_size(step,P_average,P_target);
     
      step : current step 
      P_current : Current acceptance ratio
      P_target  : preferred acceptance ratio (def=0.3);
     
      See also sippi_compute_acceptance_rate, sippi_prior_set_steplength
     

sippi\_anneal\_temperature
--------------------------

      sippi_anneal_temperature : compute annealing temperature for
      annealing type sampling
     
        %% ANNEALING (TEMPERATURE AS A FUNTION OF ITERAITON NUMBER)
         i % iteration number
     
         mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
         mcmc.anneal.i_end=100000; %  iteration number when annealing stops
         mcmc.anneal.T_begin=5; % Start temperature for annealing 
         mcmc.anneal.T_end=1; % End temperature for anneaing 
     
         mcmc.anneal.type='exp';     % Exponential temperature change
         mcmc.anneal.type='linear';  % Linear temperature change
      
      Call
        [T,mcmc]=sippi_anneal_temperature(i,mcmc);
     
      See also sippi_metropolis
     

sippi\_compute\_acceptance\_rate
--------------------------------

      sippi_compute_acceptance_rate Computes acceptance rate for the Metropolis sampler in SIPPI
     
      Call:
        P_acc=sippi_compute_acceptance_rate(acc,n_update_history);
     

sippi\_compute\_modelization\_forward\_error
--------------------------------------------

      sippi_compute_modelization_forward_error Computes an estimate of the modelization erro
     
      Computes and estimate of the Gaussian modelization error, N(dt,Ct)
      caused by the use of an imperfect forward kernel
     
      If called with only one output '[Ct]=sippi..]' then the Gaussian model is
      assumed by centered around 0, (dt{1}=0).
     
      Call
        [Ct,dt,dd]=sippi_compute_modelization_forward_error(forward_full,forward_app,prior,N,d);
     
     
      For details see:
       Hansen, T.M., Cordua, K. S., Jacobsen, B. J., and Mosegaard, K. (2014)
       Accounting for imperfect forward modeling in geophysical inverse problems - exemplified for cross hole tomography.
       Geophsyics, 79(3) H1-H21, 2014. doi:10.1190/geo2013-0215.1
     

sippi\_forward
--------------

      sippi_forward Simple forward wrapper for SIPPI
     
      Assumes that the actual forward solver has been defined by
      forward.forward_function
     
      Call:
        [d,forward,prior,data]=sippi_forward(m,forward)
     
      Optional: 
        [d,forward,prior,data]=sippi_forward(m,forward,prior)
        [d,forward,prior,data]=sippi_forward(m,forward,prior,data)
        [d,forward,prior,data]=sippi_forward(m,forward,prior,data,options)
     

sippi\_forward\_linear
----------------------

      sippi_forward_linear: 
     
      % options:
      forward.G : Linear forward operator. such that d[
                  d{id}=forward.G*m{im}(:)
                  if not set, forward.G=eye(prod(size(m{im})));
     
      forward.force_sparse [0]: Use forward.G as is (default)
                          [1]: force forward.G to be treated as a sparse matrix
     
      %% Examples
      % Define an example of a prior model
         forward.forward_function='sippi_forward_linear';
         im=1;
         prior{im}.type='FFTMA';
         prior{im}.x=[0:1:100];
         prior{im}.m0=0;
         prior{im}.Va='1 Sph(10)';
      % Define the use of the linear forward solver
         forward.forward_function='sippi_forward_linear';
     
      % Example 1: Identity, 
      %             by default an identity operator is used if the linear
      %             is not set
        m=sippi_prior(prior);
        tic;
        d=sippi_forward(m,forward);
        t1=toc;
        figure(1);
        subplot(1,2,1);plot(m{1});title('m')
        subplot(1,2,2);plot(d{1});title('d')
        suptitle(sprintf('Ex1, Identity, t=%g',t1))
     
      % Example 2: Linear Operator
        nd=prod(size(m{im}));
        forward.G=precal_cov_2d(nd,1,1,1,'1 Sph(10)');
        tic;
        d=sippi_forward(m,forward);
        t2=toc;
        figure(2);
        subplot(1,2,1);plot(m{1});title('m')
        subplot(1,2,2);plot(d{1});title('d')
        suptitle(sprintf('Ex2, t=%g',t2))
     
     
       % Example 3: Force forward.G to be sparse
        forward.force_sparse=1;
        tic;
        d=sippi_forward(m,forward);
        t3=toc;
        figure(3);
        subplot(1,2,1);plot(m{1});title('m')
        subplot(1,2,2);plot(d{1});title('d')
        suptitle(sprintf('Force sparse, t=%g',t3))
     
      See also: sippi_forward

sippi\_get\_resim\_data
-----------------------

      sippi_get_resim_data: Get conditional data for resimulation
     
      d_cond=sippi_get_resim_data(m_current,prior,ip);
      
      c_cond [n_cond,4]: col1: x, col2: y, col4: z, col4: d
     
      See also sippi_prior
     

sippi\_get\_sample
------------------

      sippi_get_sample: Get a posterior sample
     
      Call :
       [reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(working_directory,im,n_reals,skip_seq_gibbs);
     
         im: A priori model type
         n_reals: Number of realizations to return
         skip_seq_gibbs [1] Skip all realization where sequential gibbs is enabled
                        [0] Use all realization
         data: SIPPI data structure
         prior: SIPPI prior structure
         options: options structure when running sippi_metropolis
     
     
      If located in a SIPPI output folder one can simple use :
         [reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(im,n_reals);
      or
         skip_seq_gibbs=0;
         [reals,etype_mean,etype_var,reals_all,reals_ite]=sippi_get_sample(im,n_reals,skip_seq_gibbs);
     
     

sippi\_least\_squares
---------------------

      sippi_least_squares Least squares type inversion for SIPPI
     
      Call : 
         [m_reals,m_est,Cm_est]=sippi_least_squares(data,prior,forward,n_reals,lsq_type,id,im);
     
     
     
        lsq_type : 'lsq' (def), classical least squares 
                   'error_sim', simulation through error simulation
                   'visim', simulation through SGSIM of DSSIM
     

sippi\_likelihood
-----------------

      sippi_likelihood Compute likelihood given an observed dataset
     
      Call
        [logL,L,data]=sippi_likelihood(d,data);
     
     
       data{1}.d_obs [N_data,1] N_data data observations
       data{1}.d_std [N_data,1] N_data uncorrelated Gaussian STD
     
       data{1}.d_var [N_data,1] N_data uncorrelated Gaussian variances
     
     
      Gaussian modelization error, N(dt,Ct), is specified as
       data{1}.dt [N_data,1] : Bias/mean of modelization error
       data{1}.Ct [N_data,N_data] : Covariance of modelization error
     
       data{1}.Ct [1,1] : Constant Covariance of modelization error
                          imples data{1}.Ct=ones(N_data.N_data)*data{1}.Ct;
     
      data{id}.recomputeCD [default=0], if '1' then data{1}.iCD is recomputed
      each time sippi_likelihood is called. This should be used if the noise model
      changes between each call to sippi_likelihood.
     
       data{id}.full_likelihood [default=]0; if '1' the the full likelihood
       (including the determinant) is computed. This not needed if the data
       civariance is constant, but if it changes, then use
       data{id}.full_likelihood=1;
     
     
       A new type of noise model can be used as long as it is available in a 
       m file staring with 'sippi_likelihood_'. Further, it should provide the
       inputs and outputs as sippi_likelihood.m
       If a noise model has been implemented in the m-files
       sippi_likelihood_other.m
       then this can be used to evaluate the likelhood in sippi using 
       data{1}.noise_model='sippi_likelihood_other',
     
     

sippi\_mcmc\_init
-----------------

      sippi_mcmc_init Initialize McMC options for Metropolis and rejection sampling in SIPPI
     
      Call:
         options=sippi_mcmc_init(options,prior);
     

sippi\_metropolis
-----------------

      sippi_metropolis Extended Metropolis sampling in SIPPI
     
      Metropolis sampling.
        See e.g.
          Hansen, T. M., Cordua, K. S., and Mosegaard, K., 2012.
          Inverse problems with non-trivial priors - Efficient solution through Sequential Gibbs Sampling.
          Computational Geosciences. doi:10.1007/s10596-011-9271-1.
     
          Sambridge, M., 2013 - A Parallel Tempering algorithm for
          probabilistic sampling and multi-modal optimization.
     
      Call :
         [options,data,prior,forward,m_current]=sippi_metropolis(data,prior,forward,options)
      Input :
         data : sippi data structure
         prior : sippi prior structure
         forward : sippi forward structure
     
      options :
     
         options.mcmc.nite=30000;   % [1] : Number if iterations
         options.mcmc.i_sample=100; % : Number of iterations between saving model to disk
         options.mcmc.i_plot=50;  % [1]: Number of iterations between updating plots
         options.mcmc.i_save_workspace=10000;  % [1]: Number of iterations between
                                                 saving the complete workspace
         options.mcmc.i_sample=100; % : Number of iterations between saving model to disk
     
         options.mcmc.m_init : Manually chosen starting model
         options.mcmc.m_ref  : Reference known target model
     
         options.mcmc.accept_only_improvements [0] : Optimization
         options.mcmc.accept_all [0]: accepts all proposed models (ignores lilkelihood)
     
         options.txt [string] : string to be used as part of all output file names
     
         %% PERTURBATION STRATEGY
         % Perturb all model parameter all the time
         options.mcmc.pert_strategy.perturb_all=1; % Perturb all priors in each
                                                   % iteration. def =[0]
         options.mcmc.pert_strategy.perturb_all=2; % Perturb a random selection of 
                                                   % all priors in each iteration. def =[0]
     
         % Perturb one a prior type at a time, according to some frequency
         options.mcmc.pert_strategy.i_pert = [1,3]; % only perturb prior 1 and 3
         options.mcmc.pert_strategy.i_pert_freq = [2 8]; % perturb prior 3 80% of
                                                    % the time and prior 1 20%
                                                    % of the time
         % the default pertubation strategt is to select one prior model to
         % perturb at tandom for each iteration
     
     
         %% TEMPERING
         options.mcmc.n_chains=1; % set number of chains (def=1)
         options.mcmc.T=1;      % set temperature of chains [1:n_chains]
         options.mcmc.chain_frequency_jump=0.1; % probability allowing a jump
                                                %  between two chains
         %% ANNEALING (TEMPERATURE AS A FUNCTION OF ITERATION NUMBER)
         options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
         options.mcmc.anneal.i_end=100000; %  iteration number when annealing stops
         options.mcmc.anneal.T_begin=5; % Start temperature for annealing
         options.mcmc.anneal.T_end=1; % End temperature for annealing
     
         %% VERBOSITY
         The amount of text info displayed at the prompt, can be controlled by
         setenv('SIPPI_VERBOSE_LEVEL','2') % all: information on chain swapping
         setenv('SIPPI_VERBOSE_LEVEL','1') % information about seq-gibbs step update
         setenv('SIPPI_VERBOSE_LEVEL','0'); % [def] frequent update
         setenv('SIPPI_VERBOSE_LEVEL','-1'); % rare update om finish time
         setenv('SIPPI_VERBOSE_LEVEL','-2'); % indication of stop and start
         setenv('SIPPI_VERBOSE_LEVEL','-3'); % none
     
      See also sippi_rejection
     
     

sippi\_prior
------------

      sippi_prior: A priori models for SIPPI
     
      To generate a realization of the prior model defined by the prior structure use:
        [m_propose,prior]=sippi_prior(prior);
     
      To generate a realization of the prior model defined by the prior structure,
      in the vicinity of a current model (using sequential Gibbs sampling) use:
        [m_propose,prior]=sippi_prior(prior,m_current);
     
      The following types of a priori models can be used
        % two point statistics bases
        GAUSSIAN   [1D] : 1D generalized gaussian model
        UNIFORM [1D-3D] : 1D-3D uncorrelated uniform distribution
        CHOLESKY[1D-3D] : based on Cholesky decomposition
        FFTMA   [1D-3D] : based on the FFT-MA method (Multivariate Gaussian)
        VISIM   [1D-3D] : based on Sequential Gaussian and Direct Sequential simulation
        SISIM   [1D-3D] : based on Sequential indicator SIMULATION
        % multiple point based statistics
        SNESIM_STD  [1D-3D] : (SGEMS) based on a multiple point statistical model inferref from a training images. Relies in the SNESIM algorithm
        SNESIM  [1D-3D] : (GSLIB STYLE) based on a multiple point statistical model inferref from a training images. Relies in the SNESIM algorithm
     
     
     %%% SIMPLE EXAMPLE %%%
     
     % A simple 2D multivariate Gaissian based prior model based on the
     % FFT-MA method, can be defined using
        im=1;
        prior{im}.type='FFTMA';
        prior{im}.name='A SIMPLE PRIOR';
        prior{im}.x=[0:1:100];
        prior{im}.y=[0:1:100];
        prior{im}.m0=10;
        prior{im}.Va='1 Sph(10)';
        prior=sippi_prior_init(prior);
     % A realization from this prior model can be generated using
        m=sippi_prior(prior);
     % This realization can now be plotted using
        sippi_plot_prior(m,prior);
     % or
        imagesc(prior{1}.x,prior{1}.y,m{1})
     
     %%% A PRIOR MODEL WITH SEVERAL 'TYPES OF A PRIORI MODEL'
     
        im=1;
        prior{im}.type='GAUSSIAN';
        prior{im}.m0=100;
        prior{im}.std=50;
        prior{im}.norm=100;
        im=im+1;
        prior{im}.type='FFTMA';
        prior{im}.x=[0:1:100];
        prior{im}.y=[0:1:100];
        prior{im}.m0=10;
        prior{im}.Cm='1 Sph(10)';
        im=im+1;
        prior{im}.type='VISIM';
        prior{im}.x=[0:1:100];
        prior{im}.y=[0:1:100];
        prior{im}.m0=10;
        prior{im}.Cm='1 Sph(10)';
        im=im+1;
        prior{im}.type='SISIM';
        prior{im}.x=[0:1:100];
        prior{im}.y=[0:1:100];
        prior{im}.m0=10;
        prior{im}.Cm='1 Sph(10)';
        im=im+1;
        prior{im}.type='SNESIM';
        prior{im}.x=[0:1:100];
        prior{im}.y=[0:1:100];
     
        sippi_plot_prior(prior);
     
     %% Sequential Gibbs sampling
     
        All a priori model types can be perturbed, such that a new realization
        is generated in the vicinity of a current model.
        To do this Sequential Gibbs Sampling is used.
        For more information, see <a href="matlab:web('http://dx.doi.org/10.1007/s10596-011-9271-1')">Hansen, T. M., Cordua, K. S., and Mosegaard, K., 2012. Inverse problems with non-trivial priors - Efficient solution through Sequential Gibbs Sampling. Computational Geosciences</a>.
        The type of sequential Gibbs sampling can be controlled in the
        'seq_gibbs' structures, e.g. prior{1}.seq_gibbs
     
        im=1;
        prior{im}.type='SNESIM';
        prior{im}.x=[0:1:100];
        prior{im}.y=[0:1:100];
     
        [m,prior]=sippi_prior(prior);
        prior{1}.seq_gibbs.step=1; % Large step--> independant realizations
        prior{1}.seq_gibbs.step=.1; % Smaller step--> Dependant realizations
        for i=1:30;
           [m,prior]=sippi_prior(prior,m); % One iteration of Sequential Gibbs
           sippi_plot_prior(prior,m);
        end
     
      See also: sippi_prior_init, sippi_plot_prior, sippi_plot_prior_sample, sippi_prior_set_steplength.m
     
      TMH/2012
     

sippi\_prior\_cholesky
----------------------

      sippi_prior_cholesky : Cholesky type Gaussian prior for SIPPI
     
     % Example:
        ip=1;
        prior{ip}.type='cholesky';
        prior{ip}.m0=10;
        prior{ip}.Cm='.001 Nug(0) + 1 Gau(10)';
        prior{ip}.x=0:1:100;linspace(0,100,20);
        prior{ip}.y=0:1:50;linspace(0,33,30);
        [m,prior]=sippi_prior_cholesky(prior);
        sippi_plot_prior(prior,m);
     
     % Sequential Gibbs sampling
        prior{1}.seq_gibbs.step=.1;
        for i=1:100;
            [m,prior]=sippi_prior_cholesky(prior,m);
            sippi_plot_prior(prior,m);
            caxis([8 12]);drawnow;
        end
     
     % Prior covariance model
      The prior covarince model can be setup using
        prior{ip}.m0=10;
        prior{ip}.Cm='.001 Nug(0) + 1 Gau(10)';
      or
        prior{ip}.m0=10;
       and the 'Cmat' variable 'prior{ip}.Cmat' which much the contain a full
       nd X nd size covariance matrix.
       (it is computed the first the sippi_prior_cholesky is called)
     
      See also: gaussian_simulation_cholesky
     

sippi\_prior\_dsim
------------------

      sippi_prior_dsim : Direct simulation in SIPPI
     
      Example: 
       
      prior{1}.type='dsim';
      prior{1}.x=1:1:40;;
      prior{1}.y=1:1:30;;
      prior{1}.ti=channels;;
     
      m=sippi_prior(prior);
      sippi_plot_prior(prior,m);
     
     
     
      % OPTIONAL OPTIONS
     
        prior{1}.options.n_cond [int]: number of conditional points (def=5)
        prior{1}.options.n_max_ite [int]: number of maximum iterations through the TI for matching patterns (def=200)
      
        prior{1}.options.plot    [int]: [0]:none, [1]:plot cond, [2]:storing movie (def=0)
        prior{1}.options.verbose [int]: [0] no infor to screen, [1]:some info (def=1)
      
      
     
      TMH/2014
     
      See also: sippi_prior_init, sippi_prior
     

sippi\_prior\_init
------------------

      sippi_prior_init Initialize PRIOR structure for SIPPI
     
      Call
        prior=sippi_prior_init(prior);
     
      See also sippi_prior
     

sippi\_prior\_mps
-----------------

      sippi_prior_mps : prior based on MPS
     
                           Using SNESIM/ENESIM FROM 
                           https://github.com/ergosimulation/mpslib
      
     % Example:
         ip=1;
         prior{ip}.type='mps';
         prior{ip}.method='mps_snesim';
         prior{ip}.x=1:1:80;
         prior{ip}.y=1:1:80;
         prior{ip}.ti=channels;
         % prior{ip}.ti=maze;
     
         m=sippi_prior(prior);
         sippi_plot_prior(prior,m)
         figure(1);imagesc(prior{ip}.ti);axis image
     
     
     %  The specific algorithm use for MPS simulation is defined in ytje 'method' field
         prior{ip}.method='mps_snesim';         % default, same as 'mps_snesim_tree'
         prior{ip}.method='mps_snesim_tree';
         prior{ip}.method='mps_snesim_list';
         prior{ip}.method='mps_genesim';
     
         All properties for each algorithm are availale in the prior{ip}.MPS
         field
     
      See also: sippi_prior, ti, mps_cpp

sippi\_prior\_plurigaussian
---------------------------

      sippi_prior_plurigaussian: Plurigaussian type prior for SIPPI
     
     % Example:
        % PluriGaussian based on one Gaussian model / truncated Gaussian
        ip=1;
        prior{ip}.type='plurigaussian';
        prior{ip}.x=1:1:100;
        prior{ip}.y=1:1:100;
        prior{ip}.Cm='1 Gau(10)';
        prior{ip}.pg_map=[0 0 0 0 1 1 0 0 2 2 2];
     
        % PluriGaussian based on two Gaussian models 
        ip=ip+1;
        prior{ip}.type='plurigaussian';
        prior{ip}.x=1:1:100;
        prior{ip}.y=1:1:100;
        prior{ip}.pg_prior{1}.Cm=' 1 Gau(10)';
        prior{ip}.pg_prior{2}.Cm=' 1 Sph(10,90,.4)';
        prior{ip}.pg_map=[0 0 0 1 1; 1 2 0 1 1; 1 1 1 1 1];
     
        [m,prior]=sippi_prior(prior);
        sippi_plot_prior(prior,m);
     
        sippi_plot_prior_sample(prior,1,5);
        sippi_plot_prior_sample(prior,2,5);
     
     
     % Sequential Gibbs sampling
        prior{1}.seq_gibbs.step=.01;
        for i=1:100;
            [m,prior]=sippi_prior(prior,m);
            sippi_plot_prior(prior,m);
            caxis([0 2]);drawnow;
        end
     
      See also: sippi_prior, pg_transform
     

sippi\_prior\_set\_steplength
-----------------------------

      sippi_prior_set_steplength Set step length for Metropolis sampler in SIPPI
     
      Call
        prior=sippi_prior_set_steplength(prior,mcmc,im);
     

sippi\_prior\_sisim
-------------------

      sippi_prior_sisim: SISIM (SGeMS) type prior for SIPPI
     
     % Example:
         ip=1;
         prior{ip}.type='sisim';
         prior{ip}.x=1:1:80;
         prior{ip}.y=1:1:80;
         prior{ip}.Cm='1 Sph(60)';
         prior{ip}.marginal_prob=[.1 .4 .5];
         m=sippi_prior(prior);
         sippi_plot_prior(prior,m)
     
         % optionally a specific random seed can be set using
         prior{ip}.seed=1;
     
     % Sequential Gibbs sampling type 1 (box selection of pixels)
         prior{ip}.seq_gibbs.type=1;%
         prior{ip}.seq_gibbs.step=10; % resim data in 10x10 pixel grids
         [m,prior]=sippi_prior(prior);
         for i=1:10;
            [m,prior]=sippi_prior(prior,m);
            sippi_plot_prior(prior,m);
            drawnow;
         end
     
     % Sequential Gibbs sampling type 2 (random pixels)
         prior{ip}.seq_gibbs.type=2;%
         prior{ip}.seq_gibbs.step=.6; % Resim 60% of data
         [m,prior]=sippi_prior(prior);
         for i=1:10;
            [m,prior]=sippi_prior(prior,m);
            sippi_plot_prior(prior,m);
            drawnow;
         end
     
      See also: sippi_prior, sgems
     

sippi\_prior\_snesim
--------------------

      sippi_prior_snesim : SNESIM type Gaussian prior for SIPPI
     
                           Using SNESIM form
                           https://github.com/SCRFpublic/snesim-standalone
     
     % Example:
         ip=1;
         prior{ip}.type='snesim';
         prior{ip}.x=1:1:80;
         prior{ip}.y=1:1:80;
         prior{ip}.ti=channels;
         % prior{ip}.ti=maze;
     
         m=sippi_prior(prior);
         sippi_plot_prior(prior,m)
         figure(1);imagesc(prior{ip}.ti);axis image
     
     % Example: scaling and rotation
         ip=1;
         prior{ip}.type='snesim';
         prior{ip}.x=1:1:80;
         prior{ip}.y=1:1:80;
         prior{ip}.ti=channels;
         prior{ip}.scaling=[.1];
         prior{ip}.rotation=[10];
     
         m=sippi_prior(prior);
         sippi_plot_prior(prior,m)
         figure(1);imagesc(prior{ip}.ti);axis image
     
     % Hard data
        % hard data are given using either matrix of 4 columns (x,y,z,val)
        % or as a 4 column EAS file (x,y,z,val)
        d_hard=[1 1 0 0; 1 2 0 0; 2 2 0 1 ];
        prior{ip}.hard_data=d_hard;
     
        write_eas('snesim_hard.dat',d_hard);
        prior{ip}.hard_data='snesim_hard.dat';
     
     % Soft data
        % soft mush be provided as a matrix of the same size as the simulation
        % grid
        d_soft(:,:,1)=ones(80,80).*NaN
        d_soft(:,:,2)=1-d_soft(:,:,1);
        prior{ip}.soft_data_grid=d_soft;
     
     
        % Optionally the soft data can be provided as point data, in which case
        % a grid, the size of the simulation grid, with soft data values will be computed
        d_soft=[1 1 0 0.2 0.8; 1 2 0 0.1 0.9; 2 2 0  0.05 0.95];
        prior{ip}.soft_data=d_soft;
     
     % Sequential Gibbs sampling type 1 (box selection of pixels)
         prior{ip}.seq_gibbs.type=1;%
         prior{ip}.seq_gibbs.step=10; % resim data in 10x10 pixel grids
         [m,prior]=sippi_prior(prior);
         for i=1:10;
            [m,prior]=sippi_prior(prior,m);
            sippi_plot_prior(prior,m);
            drawnow;
         end
     
     % Sequential Gibbs sampling type 2 (random pixels)
         prior{ip}.seq_gibbs.type=2;%
         prior{ip}.seq_gibbs.step=.6; % Resim 60% of data
         [m,prior]=sippi_prior(prior);
         for i=1:10;
            [m,prior]=sippi_prior(prior,m);
            sippi_plot_prior(prior,m);
            drawnow;
         end
     
      See also: sippi_prior, ti
     

sippi\_prior\_snesim\_std
-------------------------

      sippi_prior_snesim_std : SNESIM_STD (SGeMS) type Gaussian prior for SIPPI
                           
                           Requires SGeMS version 2.1b, available from 
                           http://sgems.sourceforge.net/?q=node/77
      
     % Example:
         ip=1;
         prior{ip}.type='snesim_std';
         prior{ip}.x=1:1:80;
         prior{ip}.y=1:1:80;
         prior{ip}.ti=channels;
         % prior{ip}.ti=maze;
     
         m=sippi_prior(prior);
         sippi_plot_prior(prior,m)
         figure(1);imagesc(prior{ip}.ti);axis image
     
     % Hard data
        % hard data are given using either matrix of 4 columns (x,y,z,val)
        % or as a 4 column EAS file (x,y,z,val)
        d_hard=[1 1 0 0; 1 2 0 0; 2 2 0 1 ];
        prior{ip}.hard_data=d_hard;
     
        write_eas('snesim_hard.dat',d_hard);
        sgems_write_pointset('snesim_hard.sgems',d_hard);
        prior{ip}.hard_data='snesim_hard.sgems';
     
     % Example: scaling and rotation
         ip=1;
         prior{ip}.type='snesim_std';
         prior{ip}.x=1:1:80;
         prior{ip}.y=1:1:80;
         prior{ip}.ti=channels;
         prior{ip}.scaling=[.1];
         prior{ip}.rotation=[10];
     
         m=sippi_prior(prior);
         sippi_plot_prior(prior,m)
         figure(1);imagesc(prior{ip}.ti);axis image
     
     % Sequential Gibbs sampling type 1 (box selection of pixels)
         prior{ip}.seq_gibbs.type=1;%    
         prior{ip}.seq_gibbs.step=10; % resim data in 10x10 pixel grids
         [m,prior]=sippi_prior(prior);
         for i=1:10;
            [m,prior]=sippi_prior(prior,m);
            sippi_plot_prior(prior,m);
            drawnow;
         end
     
     % Sequential Gibbs sampling type 2 (random pixels)
         prior{ip}.seq_gibbs.type=2;%    
         prior{ip}.seq_gibbs.step=.6; % Resim 60% of data
         [m,prior]=sippi_prior(prior);
         for i=1:10;
            [m,prior]=sippi_prior(prior,m);
            sippi_plot_prior(prior,m);
            drawnow;
         end
     
      See also: sippi_prior, sippi_prior_snbesim, ti
     

sippi\_prior\_uniform
---------------------

      sippi_prior_uniform : Uniform prior for SIPPI
     
     
     % Example 1D uniform
        ip=1;
        prior{ip}.type='uniform';
        prior{ip}.min=10;
        prior{ip}.max=25;
        [m,prior]=sippi_prior_uniform(prior);
        sippi_plot_prior_sample(prior);
     
     % Example 10D uniform
        ip=1;
        prior{ip}.type='uniform';
        prior{ip}.x=1:1:10; % As dimensions are uncorrelated, only the lehgth
                            % of prior{ip}.x matters, not its actual values.
        prior{ip}.min=10;
        prior{ip}.max=25;
        [m,prior]=sippi_prior_uniform(prior);
        sippi_plot_prior_sample(prior);
     
     % Sequential Gibbs sampling
        prior{1}.seq_gibbs.step=.1;
        for i=1:1000;
            [m,prior]=sippi_prior(prior,m);
            mm(i)=m{1};
        end
        subplot(1,2,1);plot(mm);
        subplot(1,2,2);hist(mm);
     
      TMH/2014
     
      See also: sippi_prior_init, sippi_prior
     

sippi\_prior\_visim
-------------------

      sippi_prior_visim : VISIM type Gaussian prior for SIPPI
     
     % Example:
         ip=1;
         prior{ip}.type='visim';
         prior{ip}.x=1:1:80;
         prior{ip}.y=1:1:80;
         prior{ip}.Cm='1 Sph(60)';
         m=sippi_prior(prior);
         sippi_plot_prior(prior,m)
     
         % optionally a specific random seed can be set using
         prior{ip}.seed=1;
     
     % Sequential Gibbs sampling type 1 (box selection of pixels)
         prior{ip}.seq_gibbs.type=1;
         prior{ip}.seq_gibbs.step=10; % resim data in 10x10 pixel grids
         prior{ip}.cax=[-2 2];
         [m,prior]=sippi_prior(prior);
         for i=1:1000;
            [m,prior]=sippi_prior(prior,m);
            sippi_plot_prior(prior,m);
            drawnow;
         end
     
     % Sequential Gibbs sampling type 2 (random pixels)
         prior{ip}.seq_gibbs.type=2;%
         prior{ip}.seq_gibbs.step=.6; % Resim 60% of data
         [m,prior]=sippi_prior(prior);
         for i=1:10;
            [m,prior]=sippi_prior(prior,m);
            sippi_plot_prior(prior,m);
            drawnow;
         end
     
     % TARGET DISTRIBUTION
      clear prior
      d_target=[7 8 9 10 11 11 12];
      d_target=[7 8 9 10 14 15 20];
     
      ip=1;
      prior{ip}.type='visim';
      prior{ip}.method='sgsim';
      % prior{ip}.method='dssim';
      prior{ip}.d_target=d_target;
      prior{ip}.cax=[min(d_target) max(d_target)];
      prior{ip}.x=1:1:80;
      prior{ip}.y=1:1:80;
      prior{ip}.Cm=sprintf('%g Gau(20)',var(d_target));
      prior{ip}.Cm=sprintf('%g Gau(20)',1);
      [m,prior]=sippi_prior(prior);
      sippi_plot_prior(prior,m);
     
      prior{ip}.seq_gibbs.step=16;
      prior{ip}.seq_gibbs.type=1;
      for i=1:10;
         [m,prior]=sippi_prior(prior,m);
         sippi_plot_prior(prior,m);
         drawnow
      end
     
     
      See also: sippi_prior, visim, nscore, inscore
     

sippi\_prior\_voronoi
---------------------

      sippi_prior_voronoi: 
     
      TMH/2014
     
      Ex: 
        prior{1}.type='voronoi';
        prior{1}.x=1:1:20;
        prior{1}.y=1:1:20;
        prior{1}.cells_N_min=2;
        prior{1}.cells_N_max=100;
        prior{1}.cells_N=10;
        [m,prior]=sippi_prior(prior);
        sippi_plot_prior(prior,m);
     
      See also: sippi_prior_init, sippi_prior
     

sippi\_rejection
----------------

      sippi_rejection Rejection sampling
     
      Call :
          options=sippi_rejection(data,prior,forward,options)
     
      input arguments
     
        options.mcmc.i_plot
        options.mcmc.nite     % maximum number of iterations
        options.mcmc.logLmax [def=1]; % Maximum possible log-likelihood value
     
        options.mcmc.rejection_normalize_log = log(options.mcmc.Lmax)
     
        options.mcmc.adaptive_rejection=1, adaptive setting of maximum likelihood
                       (def=[0])
                       At each iteration logLmax will be set if log(L(m_cur)=>options.mcmc.logLmax
     
     
        options.mcmc.max_run_time_hours = 1; % maximum runtime in hours
                                             % (overrides options.mcmc.nite if needed)
     
        options.mcmc.T = 1; % Tempering temperature. T=1, implies no tempering
     
      See also sippi_metropolis
     
     

sippi\_set\_path
----------------

      sippi_set_path Set paths for running sippi
