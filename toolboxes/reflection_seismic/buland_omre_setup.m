% buland_omre_setup: setup W, A, D, matrices for linear forward operator ala Buland and Omre (2003)
%
% Call:
%    [A,D,W]=buland_omre_setup(log_vp0,log_vs0,log_rho0,ns,angle,wavelet);
%
% Input:
%    log_vp0: logarithm of background P-velocity
%    log_vs0: logarithm of background S-velocity
%    log_rho0: logarithm of background Density
%    ns: number of
%
% Output:
%    A: Linear coefficents
%    D: Differential operator
%    W: Wavelet matrix
%
function [A,D,W]=buland_omre_setup(vp0,vs0,rho0,ns,angle,wavelet);
  na=length(angle);
  
  %% A MATRIX
  vs2vp2=exp(vs0).^2./exp(vp0).^2;
  a_vp=0.5*(1+tan(pi*angle/180).^2);
  a_vs=-4*(vs2vp2)*sin(pi*angle/180).^2;
  a_rho=.5*(1+a_vs);
  Asmall=[a_vp(:) a_vs(:) a_rho(:)];

  for ia=1:na;
    A1=eye(ns).*Asmall(ia,1);
    A2=eye(ns).*Asmall(ia,2);
    A3=eye(ns).*Asmall(ia,3);

    if ia==1
      A=[A1 A2 A3];
    else
      A=[A;A1 A2 A3];
    end
  end

  %% DIFFERENTIAL MATRIX
  D_small=zeros(ns,ns);
  %for i=1:(n-1);
  %  % pad top
  %  D_small(i,[1:2]+i-1)=[-1 1];
  %end
  for i=2:ns;
    % pad bot
    D_small(i,[1:2]+i-2)=[-1 1];
  end
  D=blkdiag(D_small,D_small,D_small);

  %% WAVELET MATRIX
  if nargin>5;
    for ia=1:na;
      Wangle{ia}=zeros(ns,ns);
      if length(wavelet)==1
        Wangle{ia}=setup_wavelet_matrix(wavelet(1).data,wavelet(1).t,ns);
      else
        Wangle{ia}=setup_wavelet_matrix(wavelet(ia).data,wavelet(ia).t,ns);
      end
      % combine wavelet matrix for each angle into one big matrix
      if ia==1;
        W=Wangle{ia};
      else
        W=blkdiag(W,Wangle{ia});
      end
    end

  else
    W=[];
  end
