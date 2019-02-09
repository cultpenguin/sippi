getinunits
----------

     GETINUNITS   Get object properties in specified units
        V = GETINUNITS(H, PROP, UNITS) returns the object property
        in the specified UNITS. It will leave the 'Units' and 'FontUnits'
        property unchanged afterwards.
     
        H is the handle of the object. If it is an M-element array of handles,
        the function will return an M-by-1 cell array. PROP can be a string or
        a cell array of strings. If it is a 1-by-N or N-by-1 cell array, the
        function will return an M-by-N cell array of values. UNITS can be a
        string or a cell array. If it is a cell array, then PROP must also be a
        cell array with the same size as UNITS, and each cell element of UNITS
        corresponds to a cell element of PROP.
     
        V = GETINUNITS(H, PROP) is the same as GET(H, PROP)
     
        Examples:
          V = GETINUNITS(H, 'Position', 'Pixels')
          V = GETINUNITS(H, {'FontSize', 'Position'}, 'Normalized')
          V = GETINUNITS(H, {'FontSize', 'Position'}, {'Points', 'Pixels'})
     
        See also GET, SET

logdet
------

     LOGDET Computation of logarithm of determinant of a matrix
     
        v = logdet(A);
            computes the logarithm of determinant of A. 
     
            Here, A should be a square matrix of double or single class.
            If A is singular, it will returns -inf.
     
            Theoretically, this function should be functionally 
            equivalent to log(det(A)). However, it avoids the 
            overflow/underflow problems that are likely to 
            happen when applying det to large matrices.
     
            The key idea is based on the mathematical fact that
            the determinant of a triangular matrix equals the
            product of its diagonal elements. Hence, the matrix's
            log-determinant is equal to the sum of their logarithm
            values. By keeping all computations in log-scale, the
            problem of underflow/overflow caused by product of 
            many numbers can be effectively circumvented.
     
            The implementation is based on LU factorization.
     
        v = logdet(A, 'chol');
            If A is positive definite, you can tell the function 
            to use Cholesky factorization to accomplish the task 
            using this syntax, which is substantially more efficient
            for positive definite matrix. 
     
        Remarks
        -------
            logarithm of determinant of a matrix widely occurs in the 
            context of multivariate statistics. The log-pdf, entropy, 
            and divergence of Gaussian distribution typically comprises 
            a term in form of log-determinant. This function might be 
            useful there, especially in a high-dimensional space.       
     
            Theoretially, LU, QR can both do the job. However, LU 
            factorization is substantially faster. So, for generic
            matrix, LU factorization is adopted. 
     
            For positive definite matrices, such as covariance matrices,
            Cholesky factorization is typically more efficient. And it
            is STRONGLY RECOMMENDED that you use the chol (2nd syntax above) 
            when you are sure that you are dealing with a positive definite
            matrix.
     
        Examples
        --------
            % compute the log-determinant of a generic matrix
            A = rand(1000);
            v = logdet(A);
     
            % compute the log-determinant of a positive-definite matrix
            A = rand(1000);
            C = A * A';     % this makes C positive definite
            v = logdet(C, 'chol');
     

pg\_plot
--------

      pg_plot: plot plurigaussian transfer function
     
      See also: pg_transform
     

pg\_transform
-------------

      pg_transform: plurigaussian transformation
     
      Call: 
         [pg_d]=pg_transform(m,pg_map,pg_limits);
           m: realizations of gaussian distribtion(s)
           pg_map: Map defining pluri-Gaussian transformation 
           pg_limits: Limit values for pg_map (def:[-3 3])
     
           pd_d: Pluri-Gaussian transformed data
      
      See also: sippi_prior_plurigaussian, prior_reals_plurigaussian
     
      See also: pg_plot
     

plotboxpos
----------

     PLOTBOXPOS Returns the position of the plotted axis region
     
      pos = plotboxpos(h)
     
      This function returns the position of the plotted region of an axis,
      which may differ from the actual axis position, depending on the axis
      limits, data aspect ratio, and plot box aspect ratio.  The position is
      returned in the same units as the those used to define the axis itself.
      This function can only be used for a 2D plot.  
     
      Input variables:
     
        h:      axis handle of a 2D axis (if ommitted, current axis is used).
     
      Output variables:
     
        pos:    four-element position vector, in same units as h

quantile
--------

     QUANTILE Empirical (sample) quantiles.
     
        For vectors Q = QUANTILE(X, P) is the empirical quantiles of X for the
        probabilities in P.  The smallest observation corresponds to P = 0 and
        the largest to P = 1.  The length of Q is LENGTH(P).
     
        For matrices, QUANTILE(X, P) is the empirical quantiles of each column.
     
        In general, QUANTILE(X, P) is the empirical quantiles along the first
        non-singleton dimension.
     
        QUANTILE(X, P, DIM) returnes the quantiles along dimension DIM.
     
        This is a MATLAB version of the R `quantile' function.

setup\_Cm\_corr {#setup_Cm_corr}
---------------

      setup_Cm_corr: setup correlated covariance model structure according to 
                     Buland and Omre (2003)
     
     
      Example:
     
       mu_omre = log([3000 2250 3000]);
       var_omre = [0.0074 0.00240 .0074 ];
       cc=[1 0.8 0.7;0.8 1 0.5;0.7 0.5 1];
       Va=sprintf('0.001 Nug(0) + 0.999 Gau(20)');
       pos=[1:1:80]';
       Cm0=precal_cov(pos,pos,Va);
       [m0,Cmat]=setup_Cm_corr(Cm0,mu_omre,var_omre,cc);
     
     
      See also; sippi_prior_cholesky
     

sippi\_verbose
--------------

      sippi_verbose : displays verbose information to the console
     
      Call:
       sippi_verbose(txt,verbose)
     
      txt [string] : text to be displayed
      verbose [integer] (def=0) : increase to see more information
     
      'vlevel' must be set in the sippi_verbose.m m-file.
     
      All entries with vebose>vlevel are displayed
     
     
      entries with a higher verbose value has a higher chance of being displayed
      that entries with lower verbose values
      verbose [0] : normal (default)
              [-1] : little info 
              [-2] : no info
     
      The verbose level can be set either using a an environmental variable
      e.g. setenv('SIPPI_VERBOSE_LEVEL','1')
      or one can call sippi_verbose with a third argument, which will set the
      verbose level. To set the verbose level to 2, use:
      sippi_verbose('',0,2);
     
      The verbose level can also be set using an environmental variable:
         %% VERBOSITY
         The amount of text info displayed at the prompt, can be controlled by
         setenv('SIPPI_VERBOSE_LEVEL','2') % all: information on chain swapping
         setenv('SIPPI_VERBOSE_LEVEL','1') % information about seq-gibbs step update
         setenv('SIPPI_VERBOSE_LEVEL','0'); % [def] frequent update
         setenv('SIPPI_VERBOSE_LEVEL','-1'); % rare update om finish time
         setenv('SIPPI_VERBOSE_LEVEL','-2'); % indication of stop and start
         setenv('SIPPI_VERBOSE_LEVEL','-3'); % none
     
     
