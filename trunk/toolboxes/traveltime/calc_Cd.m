function Cd=calc_Cd(ant_pos,var_uncor,var_cor1,var_cor2,L)
% Calc_cd Setup a covariance model to account for borehole imperfections
%
% Call: Cd=calc_Cd(ant_pos,var_uncor,var_cor1,var_cor2,L)
% This function sets up a data covariance matrix that accounts for static
% (i.e. correlated) data errors.
% 
% Inputs:
% * ant_pos: A N x 4 array that contains N combinations of transmitter/source 
% and receiver positions. The first two columns are the x- and y-coordinates
% of the transmitter/source position. The last two columns are the x- and 
% y-coordiantes of the receiver position.
% * var_uncor: The variance of the uncorrelated data errors.
% * var_cor1: The variance of the correlated data errors
% related to the transmitter/source positions.
% * var_cor2: The variance of the correlated data errors
% related to the receiver positions.
% * L: The correlation length for the correlation between the individual
% transmitter/source or receiver positions using an exponential covariance 
% function. For typical static errors the correlation length is set to a 
% small number (e.g. 10^-6).
% 
% For more details and practical examples see:
% Cordua et al., 2008 in Vadose zone journal.
% Cordua et al., 2009 in Journal of applied geophysics.
%
% Knud S. Cordua (2012)

if nargin<5, L=1e-5;end
if nargin<4, var_cor2=0;end
if nargin<3, var_cor1=0;end
if nargin<2, var_uncor=0;end
if nargin<1, help(mfilename);return;end



allxss=ant_pos(:,1);
allzss=ant_pos(:,2);
allxrs=ant_pos(:,3);
allzrs=ant_pos(:,4);

% -- Diagonal elements --
Ndata=length(allxss);

C_diag=var_uncor*eye(Ndata);

% -- Off-diagonal elements --
C_offdiag1=zeros(Ndata,Ndata);
for i=1:Ndata,
   s1=sqrt((allxss(i)-allxss).^2+(allzss(i)-allzss).^2);
   s1=s1';
   C_offdiag1(i,:)=var_cor1*exp(-1*s1/L);
end   

C_offdiag2=zeros(Ndata,Ndata);
for i=1:Ndata,
   s2=sqrt((allxrs(i)-allxrs).^2+(allzrs(i)-allzrs).^2);
   s2=s2';
   C_offdiag2(i,:)=var_cor2*exp(-1*s2/L);  
end

% Sum the individual covariance matrices:
Cd=C_diag+C_offdiag1+C_offdiag2;
