function [antpos angle dist]=fwi_ant_geometry(x0,z0,x,z,dtrn,drec,maxang)

% Call: [antpos angle dist]=fwi_ant_geometry(x0,z0,x,z,dtrn,drec,maxang);
% Calculates at antennae geoemtry setup for a typical cross-borehole survey.  
% Additionally, the the angles of the rays connecting the antennae
% positions are given.
%
% * x0 is the origo x-value
% * z0 is the origo z-value
% * x is the width of the profile measured from x0
% * z is the depth of the profile measured from z0
% * dtrn is the increment between the transmitter positions
% * drec is the increment between the receiver positions 
% * maxang is the largest accepted ray angle from horisontal
% 
% * angle is the angles of the positions listed in antpos (in degrees)
% * dist is the distance between the transmitter and receiver position
%
% All units are in meter
%
% Knud S. Cordua (2008)


% From degrees to radians:
maxang=maxang*pi/180;

Ntrn=length(z0:dtrn:(z+z0));
Nrec=length(z0:drec:(z+z0));
z_trn=z0:dtrn:(z+z0);
z_rec=z0:drec:(z+z0);

% For left located transmitter positions
x_trn=x0;
x_rec=x0+x;
c=0;
for i=1:Ntrn
    for j=1:Nrec
        ang=atan(abs(z_trn(i)-z_rec(j))/x);
        if abs(ang)<=maxang
            c=c+1;
            angle(c)=ang;
            antpos(c,1)=x_trn;
            antpos(c,2)=z_trn(i);
            antpos(c,3)=x_rec;
            antpos(c,4)=z_rec(j);
            dist(c)=sqrt((z_trn(i)-z_rec(j))^2+(x_trn-x_rec)^2);
        end
    end
end
    
% For right located transmitter positions
x_trn=x0+x;
x_rec=x0;
for i=1:Ntrn
    for j=1:Nrec
        ang=atan(abs(z_trn(i)-z_rec(j))/x);
        if abs(ang)<=maxang
            c=c+1;
            angle(c)=ang;
            antpos(c,1)=x_trn;
            antpos(c,2)=z_trn(i);
            antpos(c,3)=x_rec;
            antpos(c,4)=z_rec(j);
            dist(c)=sqrt((z_trn(i)-z_rec(j))^2+(x_trn-x_rec)^2);
        end
    end
end



