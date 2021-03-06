sippi\_colormap
---------------

      sippi_colormap Default colormap for sippi
     
      Call :
        sippi_colormap; % the same as sippi_colormap(3);
     
      or :
        sippi_colormap(1) - Red Green Black
        sippi_colormap(2) - Red Green Blue Black
        sippi_colormap(3) - Jet
        sippi_colormap(4) - Parula
        sippi_colormap(5) - Geosoft
     

sippi\_get\_posterior\_data
---------------------------

      sippi_get_posterior_data: load all data stored in mat-file
     
      Call: 
     
       [data,prior,options,mcmc]=sippi_get_posterior_data(folder_name);
       [data,prior,options,mcmc]=sippi_get_posterior_data(output_stuct);
     

sippi\_plot\_current\_model
---------------------------

      sippi_plot_current_model Plots the current model during Metropolis sampling
     
      Call :
        sippi_plot_current_model(mcmc,data,d,m_current,prior,options);
     

sippi\_plot\_data
-----------------

      sippi_plot_data: Plot data response
     
      Call:
         sippi_plot_data(d,data);
     
      sippi_plot_data provides a very simple way to plot data. 
      A more appropriate data plot can be implemented by implementing a new
      mfile called "sippi_plot_data" and add it the Matlab path before the 
      main SIPPI folders
     
      A maximum of options.plot.plot_data_max_data (def:=5) data is plotted.
      Change the value in sippi_plot_defaults.m
     
      A specific m-file for handling plotting for a specfic type of data can
      be implemented, and can be called instead of sippi_plot_data, if the name
      of the m-file is set in the options.sippi_plot_data_functions is set, as
      e.g:
         options.sippi_plot_data_function='sippi_plot_data_gpr';
      then 
         sippi_plot_data(d,data,[],options)
         sippi_plot_data(d,data,1:length(d),options); 
      will be equivalent to call
         sippi_plot_data_gpr(d,data); 
     
      
     
      
      See also sippi_plot_defaults, sippi_plot_data_gpr
     

sippi\_plot\_defaults
---------------------

      sippi_plot_defaults: Sets default options for ploting (such as fontsize)
     
      Call :
        options==sippi_plot_defaults(options);
     
        % ALWAYS USE DEFULT SETTING (overrules options.axis)
        overrule=1; % {default overrule=0)
        options==sippi_plot_defaults(options,overrule);
     
      See also: sippi_plot_posterior, sippi_plot_posterior_2d_marg
     

sippi\_plot\_loglikelihood
--------------------------

      sippi_plot_loglikelihood Plot loglikelihood time series
     
      Call : 
         acc=sippi_plot_loglikelihood(logL,i_acc,N,itext)
     

sippi\_plot\_model
------------------

sippi\_plot\_movie
------------------

      sippi_plot_movie plot movie of prior and posterior realizations
     
      Call :
        sippi_plot_movie(fname);
        sippi_plot_movie(fname,im_array,n_frames,skip_burnin,i_chain);
           fname : name of folder with results (e.g. options.txt)
           im_array : array of indexes of model parameters to make into movies
           n_frames [200] : number of frames in movie
           skip_burnin [200] : start movie after burn_in;
           i_chain[1]: make movie of chain number 'i_chain' (new 22/05/2014)
     
      Ex:
      sippi_plot_movie('20130812_Metropolis');
      sippi_plot_movie(options.txt);
     
      %% 1000 realization including burn-in, for prior number 1
      sippi_plot_movie('20130812_Metropolis',1,1000,0);
     
      Using options.plot.skip_seq_gibbs=1, (set in sippi_plot_defaults)
      removes realizations obtained using sequential Gibbs sampling 
      (equivalent to setting skip_burnin=1)
     
      See also: sippi_plot_defaults
     

sippi\_plot\_posterior
----------------------

      sippi_plot_posterior Plot statistics from posterior sample
     
      Call :
         sippi_plot_posterior(fname,im_arr,prior,options,n_reals);
     
      See also sippi_plot_prior
     

sippi\_plot\_posterior\_2d\_marg
--------------------------------

      sippi_plot_posterior_2d_marg: plots 2D posterior marginal distributions
     
      Call:
         [options,reals_all]=sippi_plot_posterior_2d_marg(options,prior,data,fname);
     
      See also: sippi_plot_posterior
     

sippi\_plot\_posterior\_data
----------------------------

      sippi_plot_posterior_data: plots posterior data and noise relaizations
     
      Call 
         [options]=sippi_plot_posterior_data(options,prior,data,forward);
            
      See also: sippi_plot_posterior
     

sippi\_plot\_posterior\_loglikelihood
-------------------------------------

      sippi_plot_posterior_loglikelihod : plots log(L) and autorreation of log(L)
     
      Call:
          sippi_plot_posterior_loglikelihood; % when located in an output folder
                                              % generated by SIPPI
     
          sippi_plot_posterior_loglikelihood(foldername); % Where 'foldername'
                                              % is a folder generated by SIPPI
     
     
          sippi_plot_posterior_loglikelihood(options); % where options is the
                             % output of sippi_rejection or sippi_metropolis
     
     
          options=sippi_plot_posterior_loglikelihood(options,prior,data,mcmc,fname);
     
     
      See also: sippi_plot_posterior
     

sippi\_plot\_posterior\_sample
------------------------------

      sippi_plot_posterior_sample: plots posterior sample statistics
     
      Call
         [options]=sippi_plot_posterior_sample(options,prior,data,forward);
     
      See also: sippi_plot_posterior
     

sippi\_plot\_prior
------------------

      sippi_plot_prior Plot a 'model', i.e. a realization of the prior model
     
     
      Call :
        sippi_plot_prior(prior,m,im_array);
     
        prior : Matlab structure for SIPPI prior model
        m : Matlab structure for SIPPI realization
        im_array : integer array of type of models to plot (typically 1)
     
     
      Example
        m=sippi_prior(prior);
        sippi_plot_prior(prior,m);
     
        m=sippi_prior(prior);
        sippi_plot_prior(prior,m,2);
     
      See also sippi_plot_prior, sippi_prior

sippi\_plot\_prior\_movie
-------------------------

      sippi_plot_prior_movie: creates a movie file a random walk in the prior
      
      Call:
         sippi_plot_prior_movie(prior,n_frames,options,im_array);
     
      See also: sippi_plot_movie
     

sippi\_plot\_prior\_sample
--------------------------

      sippi_plot_prior_sample Plot a sample of the prior in SIPPI
     
      Call :
         sippi_plot_prior_sample(prior,im_array,n_reals,cax,options);
     
       See also sippi_plot_posterior, sippi_plot_prior, sippi_prior
     

sippi\_plot\_set\_axis
----------------------

      sippi_plot_set_axis
      see also sippi_plot_defaults

wiggle
------

      wiggle : plot wiggle/VA/image plot
     
      Call
         wiggle(Data); % wiggle plot
         wiggle(Data,scale); % scaled wiggle plot
         wiggle(x,t,Data); % wiggle plt
         wiggle(x,t,Data,'VA') % variable Area (pos->black;neg->transp)
         wiggle(x,t,Data,'VA2') % variable Area (pos->black;neg->red)
         wiggle(x,t,Data,'wiggle',scale); % Scaled wiggle
         wiggle(x,t,Data,'wiggle',scale,showmax); % Scaled wiggle and max
                                                    showmax traces.
         wiggle(x,t,Data,'wiggle',scale,showmax,plimage); % wiggle + image
         wiggle(x,t,Data,'wiggle',scale,showmax,plimage,caxis); % wiggle +
                                                                  scaled image
     
      Data : [nt,ntraces]
      x : [1:ntraces] X axis (ex [SegyTraceheaders.offset])
      t : [1:nt] Y axis
      style : ['VA'] : Variable Area
              ['wiggle'] : Wiggle plot
      scale : scaling factor, can be left empty as []
      showmax [scalar] : max number of traces to show on display [def=100]
      plimage [0/1] : Show image beneath wiggles [def=0];
      caxis [min max]/[scalar] : amplitude range for colorscale
     
