% correlated_traveltime_tomography_errors
%
% Call
%    [Ct, Cs, Cr, C]=correlated_traveltime_tomography_errors(sources,receivers,cTrans,cRec,cUnc)
%
%    sources [nd,ndim]
%    receivers [nd,ndim]
%    cTrans: variance of transmitter static errors
%    cRec: variance of Receiver static errors
%    cUnc: variance of uncorrelated errors
% 
% Based on 
% Cordua, Knud S., Majken C. Looms, and Lars Nielsen. "Accounting for correlated data errors during inversion of cross-borehole ground penetrating radar data." Vadose Zone Journal 7.1 (2008): 263-271.
% doi: https://doi.org/10.2136/vzj2007.0008
%
% TMH/2024
% 
function [Ct, Cs, Cr, C]=correlated_traveltime_tomography_errors(sources,receivers,cTrans,cRec,cUnc)
if nargin<2, 
    help(mfilename);
    return
end
if nargin<3, cTrans = 1; end
if nargin<4, cRec = 0; end
if nargin<5, cUnc = 0; end


nd=size(sources,1);
C=eye(nd).*cUnc;
Cs=zeros(nd,nd);
Cr=zeros(nd,nd);

%
for id1=1:nd
for id2=1:nd
    dS=edist(sources(id1,:)-sources(id2,:));
    dR=edist(receivers(id1,:)-receivers(id2,:));
    if dS==0; Cs(id1,id2)=Cs(id1,id2)+cTrans;end
    if dR==0; Cr(id1,id2)=Cr(id1,id2)+cRec;end

end
end

% Combine sources of errors
Ct=C+Cs+Cr;
