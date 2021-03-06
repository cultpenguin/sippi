sippi\_forward\_covariance\_inference
-------------------------------------

      sippi_forward_covariance_inference : Probabilitsic covariance inference
     
      Call :
       [d,forward,prior,data]=sippi_forward_covariance_inference(m,forward,prior,data,id,im)
     
     
       forward.pos_known : [x' y' z'],  [ndata,ndim] with position of observed data
       forward.G : Forward operator
     
       Prior covariance model, N(m0,Cm) is chosen as
       forward.m0 : initial mean model
       forward.Cm/forward.Va : inital covariance model
     
       prior{im}.m0
       prior{im}.type (round(type) is itype)
       prior{im}.range_1
       prior{im}.range_2
       prior{im}.range_3
       prior{im}.ang_1
       prior{im}.ang_2
       prior{im}.ang_3
       prior{im}.sill
       prior{im}.nugget_fraction
     
      See Hansen et al. (2014). A general probabilistic approach for inference of Gaussian model parameters from noisy data of point and volume support. 
      Mathematical Geosciences 47(7), pp 843-865. published online 09-2014. doi:10.1007/s11004-014-9567-5 
     
      See also: sippi_forward
     
